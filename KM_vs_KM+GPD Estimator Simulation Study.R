## Load the libraries

library(survival)
library(parallel)
library(foreach)
library(doParallel)

# Set up parallel processing
num_cores <- detectCores() - 1  # Leave one core free
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Export necessary packages to workers
clusterEvalQ(cl, {
  library(survival)
})

cat("Using", num_cores, "cores for parallel processing\n\n")

# Calculate Tm and M for censored data simulation
calc_censoring_params <- function(A, B, k = 1, theta = 3) {
  S_T <- function(t) 1 - pgamma(t, shape = k, scale = theta)
  
  if (A == 0 && B == 0) {
    return(list(Tm = Inf, M = Inf))
  }
  
  if (A == 0) {
    Tm <- Inf
    M <- integrate(S_T, 0, Inf)$value / B
    return(list(Tm = Tm, M = M))
  }
  
  if (B == 0) {
    Tm <- qgamma(1 - A, shape = k, scale = theta)
    return(list(Tm = Tm, M = Inf))
  }
  
  # Both types - solve for Tm
  f_Tm <- function(Tm) {
    int_val <- integrate(S_T, 0, Tm)$value
    S_T(Tm) * (int_val / B - Tm) - A * int_val / B
  }
  
  f_prime <- function(Tm) {
    h <- 1e-5
    (f_Tm(Tm + h) - f_Tm(Tm - h)) / (2 * h)
  }
  
  Tm <- qgamma(1 - A, shape = k, scale = theta)
  for (i in 1:100) {
    Tm <- Tm - f_Tm(Tm) / f_prime(Tm)
    if (abs(f_Tm(Tm)) < 1e-8) break
  }
  
  M <- integrate(S_T, 0, Tm)$value / B
  return(list(Tm = Tm, M = M))
}

# Generate censored data
generate_data <- function(n, A, B, k = 1, theta = 3) {
  true_times <- rgamma(n, shape = k, scale = theta)
  
  # Handle no censoring case
  if (A == 0 && B == 0) {
    return(data.frame(time = true_times, event = rep(1, n)))
  }
  
  params <- calc_censoring_params(A, B, k, theta)
  
  # Handle only administrative censoring
  if (B == 0) {
    obs_times <- pmin(true_times, params$Tm)
    status <- as.numeric(true_times <= params$Tm)
    return(data.frame(time = obs_times, event = status))
  }
  
  # Handle random censoring (with or without administrative)
  censor_times <- runif(n, 0, params$M)
  
  if (A == 0) {
    # Only random censoring
    obs_times <- pmin(true_times, censor_times)
    status <- as.numeric(true_times <= censor_times)
  } else {
    # Both types of censoring
    obs_times <- pmin(true_times, censor_times, params$Tm)
    status <- as.numeric(true_times <= censor_times & true_times <= params$Tm)
  }
  
  return(data.frame(time = obs_times, event = status))
}

# GPD censored log-likelihood
gpd_loglik <- function(params, t, delta) {
  xi <- params[1]
  sigma <- params[2]
  
  if (sigma <= 0) return(-Inf)
  
  z <- 1 + xi * t / sigma
  if (any(z <= 0)) return(-Inf)
  
  if (abs(xi) < 1e-8) {
    logf <- -log(sigma) - t / sigma
    logs <- -t / sigma
  } else {
    logf <- -log(sigma) - (1/xi + 1) * log(z)
    logs <- -(1/xi) * log(z)
  }
  
  return(sum(delta * logf + (1 - delta) * logs))
}


# Fitting Hybrid Model (KM till threshold and GPD beyond that)

fit_hybrid <- function(data, threshold_pct = 0.8) {
  
  km <- survfit(Surv(time, event) ~ 1, data = data)
  
  # Extract event times
  event_times <- data$time[data$event == 1]
  if (length(event_times) < 5) return(NULL)  # Too few events to fit GPD
  
  # Maximum observed event time
  Tm <- max(event_times)
  
  # Threshold for GPD fitting
  threshold <- quantile(event_times, probs = threshold_pct)
  
  # Survival probability at Tm
  S_Tm <- min(km$surv[km$time <= Tm])
  
  # Tail data above threshold for GPD fitting
  tail_data <- data[data$time > threshold, ]
  if (nrow(tail_data) < 3) return(NULL)  # Not enough tail points
  tail_data$time <- tail_data$time - threshold  # shift for GPD fitting
  
  # Fit GPD using MLE
  start <- c(0.1, sd(tail_data$time))
  opt <- optim(
    par = start,
    fn = function(p) -gpd_loglik(p, tail_data$time, tail_data$event),
    method = "L-BFGS-B",
    lower = c(-Inf, 1e-6),
    upper = c(Inf, Inf)
  )
  if (opt$convergence != 0) return(NULL)
  
  # Return fitted components
  return(list(
    km = km,
    threshold = threshold,
    Tm = Tm,
    S_Tm = S_Tm,
    xi = opt$par[1],
    sigma = opt$par[2]
  ))
}


# KM Mean: Calculate restricted mean survival from KM curve

calc_km_mean <- function(data) {
  km <- survfit(Surv(time, event) ~ 1, data = data)
  
  # Include time 0 with survival probability 1
  times <- c(0, km$time)
  survs <- c(1, km$surv)
  
  # Integrate stepwise (area under KM curve)
  area <- 0
  for (i in 1:(length(times) - 1)) {
    area <- area + (times[i + 1] - times[i]) * survs[i]
  }
  
  return(area)
  
  
  # Hybrid Mean: KM up to Tm (Max Observed time) + GPD tail beyond Tm
  
  calc_hybrid_mean <- function(data, threshold_pct = 0.8) {
    
    # Check if there is censoring
    if (all(data$event == 1)) {
      return(calc_km_mean(data))  # No censoring, KM mean is exact
    }
    
    # Otherwise, fit hybrid
    fit <- fit_hybrid(data, threshold_pct)
    if (is.null(fit)) {
      return(calc_km_mean(data))
    }
    
    # KM contribution
    km_times <- fit$km$time[fit$km$time <= fit$Tm]
    km_survs <- fit$km$surv[fit$km$time <= fit$Tm]
    times <- c(0, km_times, fit$Tm)
    survs <- c(1, km_survs)
    km_area <- sum(diff(times) * survs)
    
    # Tail contribution
    if (fit$xi >= 1) return(NA)  # Mean doesn't exist
    sigma_Tm <- fit$sigma + fit$xi * (fit$Tm - fit$threshold)
    gpd_mean_tail <- sigma_Tm / (1 - fit$xi)
    
    total_mean <- km_area + fit$S_Tm * gpd_mean_tail
    return(total_mean)
  }
  
  
  # Single simulation function
  one_sim <- function(n, A, B, k = 1, theta = 3) {
    data <- generate_data(n, A, B, k, theta)
    true_mean <- k * theta
    
    km_est <- tryCatch(calc_km_mean(data), error = function(e) NA)
    hybrid_est <- tryCatch(calc_hybrid_mean(data), error = function(e) NA)
    
    km_error <- km_est - true_mean
    hybrid_error <- hybrid_est - true_mean
    
    return(c(km_est = km_est, hybrid_est = hybrid_est,
             km_error = km_error, hybrid_error = hybrid_error))
  }
  
  # Running parallel simulations
  n_sims <- 200
  n_sample <- 500
  k <- 1
  theta <- 3
  true_mean <- k * theta
  
  scenarios <- data.frame(
    name = c("No censoring", "Light (20%)", "Moderate (30%)", "Heavy (50%)"),
    A = c(0.00, 0.10, 0.15, 0.25),
    B = c(0.00, 0.10, 0.15, 0.25)
  )
  
  cat("Running parallel simulations...\n")
  cat("Sample size:", n_sample, "\n")
  cat("Replications:", n_sims, "\n")
  cat("True mean:", true_mean, "\n\n")
  
  # Creating output directory
  if (!dir.exists("simulation_results")) {
    dir.create("simulation_results")
  }
  
  # Store all raw results
  all_raw_results <- list()
  
  for (i in 1:nrow(scenarios)) {
    cat("Scenario", i, ":", scenarios$name[i], "\n")
    
    start_time <- Sys.time()
    
    # Parallel loop - export all needed functions and data
    results <- foreach(j = 1:n_sims, .combine = rbind,
                       .packages = c("survival"),
                       .export = c("generate_data", "calc_censoring_params",
                                   "gpd_loglik", "fit_hybrid", 
                                   "calc_km_mean", "calc_hybrid_mean",
                                   "one_sim", "n_sample", "k", "theta")) %dopar% {
                                     one_sim(n_sample, scenarios$A[i], scenarios$B[i], k, theta)
                                   }
    
    end_time <- Sys.time()
    time_taken <- end_time - start_time
    
    cat("  Time taken:", round(time_taken, 2), attr(time_taken, "units"), "\n")
    
    # Store raw results
    results_df <- data.frame(
      scenario = scenarios$name[i],
      scenario_num = i,
      A = scenarios$A[i],
      B = scenarios$B[i],
      sim_num = 1:n_sims,
      km_est = results[, 1],
      hybrid_est = results[, 2],
      km_error = results[, 3],
      hybrid_error = results[, 4]
    )
    
    all_raw_results[[i]] <- results_df
    
    # Save individual scenario CSV
    filename <- paste0("simulation_results/scenario_", i, "_", 
                       gsub(" ", "_", scenarios$name[i]), ".csv")
    write.csv(results_df, filename, row.names = FALSE)
    cat("  Saved to:", filename, "\n\n")
  }
  
  # Combine all results
  combined_results <- do.call(rbind, all_raw_results)
  write.csv(combined_results, "simulation_results/all_simulations.csv", 
            row.names = FALSE)
  
  # Summary statistics
  cat("\n=== SUMMARY RESULTS ===\n")
  cat("True mean survival time:", true_mean, "\n\n")
  
  summary_table <- data.frame(
    Scenario = character(),
    A = numeric(),
    B = numeric(),
    KM_Mean = numeric(),
    Hybrid_Mean = numeric(),
    KM_MAE = numeric(),
    Hybrid_MAE = numeric(),
    KM_SD = numeric(),
    Hybrid_SD = numeric(),
    KM_Min = numeric(),
    KM_Max = numeric(),
    Hybrid_Min = numeric(),
    Hybrid_Max = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:nrow(scenarios)) {
    scenario_data <- all_raw_results[[i]]
    
    summary_table[i, ] <- list(
      Scenario = scenarios$name[i],
      A = scenarios$A[i],
      B = scenarios$B[i],
      KM_Mean = mean(scenario_data$km_est, na.rm = TRUE),
      Hybrid_Mean = mean(scenario_data$hybrid_est, na.rm = TRUE),
      KM_MAE = mean(abs(scenario_data$km_error), na.rm = TRUE),
      Hybrid_MAE = mean(abs(scenario_data$hybrid_error), na.rm = TRUE),
      KM_SD = sd(scenario_data$km_error, na.rm = TRUE),
      Hybrid_SD = sd(scenario_data$hybrid_error, na.rm = TRUE),
      KM_Min = min(scenario_data$km_est, na.rm = TRUE),
      KM_Max = max(scenario_data$km_est, na.rm = TRUE),
      Hybrid_Min = min(scenario_data$hybrid_est, na.rm = TRUE),
      Hybrid_Max = max(scenario_data$hybrid_est, na.rm = TRUE)
    )
  }
  
  print(summary_table)
  
# Save summary table
write.csv(summary_table, "simulation_results/summary_statistics.csv", 
            row.names = FALSE)
  
cat("\nAll results saved to 'simulation_results/' directory\n")
  
# Stop cluster
stopCluster(cl)
cat("\nParallel processing completed.\n")