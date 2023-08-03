fxn_glm_vanilla <- function(Y, X_orig, X_aug) {

  fitted_model <- glm(Y ~ X_orig + X_aug,
                      family = "binomial");
  conv <- (max(summary(fitted_model)$coefficients[,2]) < 1e3);
  
  beta_hat <- coef(fitted_model)[-1]
  names(beta_hat) <- 
    c(paste0("X_orig", 1:ncol(X_orig)), paste0("X_aug", 1:1:ncol(X_aug)))
  
  return(list(
    intercept_hat = coef(fitted_model)[1],
    beta_hat = beta_hat,
    iter = NA,
    converge = conv, 
    singular_hessian = NA,
    final_diff = NA,
    objective_function = c("total" = NA)))
  
}