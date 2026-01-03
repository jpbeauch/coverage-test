source("R/simulate_model.R")

#' Adjust theta parameter to make variances add up to 1.
#'
#' @param theta Coefficient vector for covariates.
#' @param var_v Variance component for v.
#' @param var_epsilon Variance component for epsilon.
#' @param vcov_w Variance-covariance matrix for w.
#' @return Adjusted theta vector scaled to target variance.
adjust_theta <- function(theta, var_v, var_epsilon, vcov_w) {
  var_theta_w <- as.vector(t(theta) %*% vcov_w %*% theta)
  target_var <- 1 - var_epsilon - var_v
  theta <- theta * sqrt(target_var / var_theta_w)
  return(theta)
}

#' Run bootstrap simulations for hapr model.
#'
#' @param n_bootstraps Number of bootstrap iterations.
#' @param n_observations Number of observations per simulation.
#' @param beta_g Coefficient for genetic factor.
#' @param beta_w Coefficient vector for covariates.
#' @param lp_mean mean of linear predictor used to set beta_constant.
#' @param theta Adjusted theta parameter vector.
#' @param var_v Variance component for v.
#' @param var_epsilon Variance component for epsilon.
#' @param e_w Expected value vector for w.
#' @param vcov_w Variance-covariance matrix for w.
#' @param model_type Model type ("lm", "probit", or "cox").
#' @return List with betas_df and se_df data frames containing bootstrap estimates.
run_bootstrap <- function(n_bootstraps, n_observations, beta_g, beta_w, theta, var_v, var_epsilon, e_w, vcov_w, model_type, beta_intercept = 0, cox_censoring_time = 10, cox_median_event_probability = 0.5) {
  improvement_ratio <- 1 / (1 - var_epsilon)
  
  beta_estimates_list <- list()
  se_estimates_list <- list()
  covariate_names <- names(beta_w)
    
  for (i in 1:n_bootstraps) {
    # Simulate dataset based on model type
    simulated_dataset <- simulate_model(n_observations, beta_g, beta_w, theta, var_v, var_epsilon, e_w, vcov_w, model_type, beta_intercept, cox_censoring_time, cox_median_event_probability)
    
    # Fit model
    fit <- hapr::hapr(
      y = simulated_dataset$y,
      gc = simulated_dataset$gc,
      w = simulated_dataset$w,
      model_type = model_type,
      improvement_ratio = improvement_ratio
    )
    
    # Extract estimates
    beta_hat <- fit$coefficients$beta
    se_hat <- fit$standard_errors
    
    # Store as list elements (each is a named vector)
    beta_estimates_list[[i]] <- beta_hat
    se_estimates_list[[i]] <- se_hat
  }
  
  # Convert to data frames - each list element becomes a row
  betas_df <- do.call(rbind, lapply(beta_estimates_list, function(x) as.data.frame(t(x))))
  betas_df$bootstrap_id <- 1:n_bootstraps
  
  se_df <- do.call(rbind, lapply(se_estimates_list, function(x) as.data.frame(t(x))))
  se_df$bootstrap_id <- 1:n_bootstraps
  
  list(
    betas_df = betas_df,
    se_df = se_df
  )
}


#########################
#' Create coverage table from bootstrap results.
#'
#' @param betas_df Data frame with bootstrap beta estimates (one row per bootstrap).
#' @param se_df Data frame with bootstrap standard errors (one row per bootstrap).
#' @param true_beta_g True value for genetic factor (gf) coefficient.
#' @param true_beta_w Named vector of true values for covariate coefficients.
#' @return Data frame with columns: variable, true_beta, avg_beta_hat, stdev_beta_hat, mean_se, coverage_pct.
coverage_table <- function(betas_df, se_df, true_beta_g, true_beta_w) {
  confidence_level <- 0.95
  z_multiplier <- qnorm(1 - (1 - confidence_level) / 2)
  
  # Helper function to calculate coverage for a single variable
  calc_var_stats <- function(var_name, true_val) {
    beta_hats <- betas_df[[var_name]]
    se_hats <- se_df[[var_name]]
    
    # Remove bootstrap_id column if present
    if (var_name == "bootstrap_id") {
      return(NULL)
    }
    
    # Calculate statistics
    avg_beta_hat <- mean(beta_hats, na.rm = TRUE)
    stdev_beta_hat <- sd(beta_hats, na.rm = TRUE)
    mean_se <- mean(se_hats, na.rm = TRUE)
    
    # Calculate coverage (percentage of times true value is in 95% CI)
    ci_lower <- beta_hats - z_multiplier * se_hats
    ci_upper <- beta_hats + z_multiplier * se_hats
    inside_ci <- sum(true_val >= ci_lower & true_val <= ci_upper, na.rm = TRUE)
    n_valid <- sum(!is.na(beta_hats) & !is.na(se_hats))
    coverage_pct <- if (n_valid > 0) 100 * inside_ci / n_valid else NA
    
    data.frame(
      variable = var_name,
      true_beta = true_val,
      avg_beta_hat = avg_beta_hat,
      stdev_beta_hat = stdev_beta_hat,
      mean_se = mean_se,
      coverage_pct = coverage_pct,
      stringsAsFactors = FALSE
    )
  }
  
  # Calculate stats for genetic factor (gf)
  results_list <- list()
  
  # Add gf if it exists in betas_df
  if ("gf" %in% colnames(betas_df)) {
    results_list[["gf"]] <- calc_var_stats("gf", true_beta_g)
  }
  
  # Add each variable from true_beta_w
  for (var_name in names(true_beta_w)) {
    if (var_name %in% colnames(betas_df)) {
      results_list[[var_name]] <- calc_var_stats(var_name, true_beta_w[[var_name]])
    }
  }
  
  # Combine into single data frame
  do.call(rbind, results_list)
}

#########################