#Pedro Orozcod del Pino+
#Data Integration - Exploiting Ratios of Parameter Estimates from a Reduced External Model
# 

#'input:
#'
#'Y is the response variable (for the current study)
#'p is the number of original covariates
#'q is the number of added covariates
#'n_curr is the sample size of the current study
#'X_orig is the original covariates in the current study
#'X_aug the added covariates in the current study
#'theta_tilde_no_intercept effect size estimates from the historical study for the original covariates (drop the intercept)



ratios <- function(Y, X_orig, X_aug, theta_tilde_no_intercept) {
  
  p <- ncol(X_orig)
  q <- ncol(X_aug)
  dimnames(X_orig) <-
    dimnames(X_aug) <- NULL
  
  #' Step 1: calculate residuals
  mod <- lm(X_aug ~ X_orig)
  gamma <- t(coef(mod))
  if(p == 1) {
    gamma_names = "X_orig"
  } else {
    gamma_names <- paste0("X_orig", 1:p)
  }
  
  W <- resid(mod)
  
  #' Step 2: fit the glm model 
  theta.X <- X_orig %*% theta_tilde_no_intercept
  fit <- glm(Y ~ theta.X + W, family = "binomial")
  alpha <- fit$coefficients
  if(q == 1) {
    alphaW_names = "W"
  } else {
    alphaW_names <- paste0("W", 1:q)
  }
  
  conv = (max(summary(fit)$coefficients[,"Std. Error"]) < 1e3);
  
  #' Step 3: obtain final estimates
  
  intercept_hat <- alpha["(Intercept)"] - sum(alpha[alphaW_names] * gamma[,"(Intercept)"]) #intercept estimate
  
  beta_hat <- 
    c(alpha["theta.X"] * theta_tilde_no_intercept - 
        alpha[alphaW_names] %*% gamma[,gamma_names],
      alpha[alphaW_names])
  names(beta_hat) <- 
    c(paste0("X_orig", 1:p), paste0("X_aug", 1:q))
  
  
  return(list(intercept_hat = intercept_hat,
              beta_hat = beta_hat,
              converge = conv, 
              iter = NA,
              singular_hessian = NA,
              final_diff = NA, 
              objective_function = c("total" = NA)))
  
  
}
