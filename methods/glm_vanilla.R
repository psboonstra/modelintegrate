fxn_glm_vanilla <- function(Y, X) {
  
  fitted_model <- glm(Y ~ ., data = data.frame(Y, X), family = "binomial");
  conv <- (max(summary(fitted_model)$coefficients[,2]) < 1e3);
  
  beta_hat <- coef(fitted_model)[-1]
  
  return(list(
    intercept_hat = coef(fitted_model)[1],
    beta_hat = beta_hat,
    converge = conv, 
    message = NA, 
    iter = NA,
    singular_hessian = NA,
    final_diff = NA,
    objective_function = c("total" = NA)))
  
}