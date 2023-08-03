fxn_glm_bayes <- function(Y, X_orig, X_aug, seed = sample.int(.Machine$integer.max, 1)) {
  
  temp_data <- data.frame(Y, X_orig, X_aug)
  
  fitted_model <- stan_glm(Y ~ .,
                           data = temp_data,
                           family = "binomial", 
                           prior = student_t(df = 3),
                           warmup = 1e3, 
                           iter = 2e3,
                           chains = 2,
                           cores = parallel::detectCores(),
                           verbose = FALSE, 
                           seed = seed, 
                           refresh = 0);
  
  beta_hat <- coef(fitted_model)[-1]
  names(beta_hat) <- 
    c(paste0("X_orig", 1:ncol(X_orig)), paste0("X_aug", 1:1:ncol(X_aug)))
  
  return(list(
    intercept_hat = coef(fitted_model)["(Intercept)"],
    beta_hat = beta_hat,
    iter = NA,
    converge = NA, 
    singular_hessian = NA,
    final_diff = NA,
    objective_function = c("total" = NA)))
  
}