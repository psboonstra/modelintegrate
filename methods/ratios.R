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
#'theta_tilde effect size estimates from the historical study for the original covariates (drop the intercept)

ratios <- function(Y, X, theta_tilde) {
  
  p <- length(theta_tilde)
  q <- ncol(X) - p
  orig_names <- names(theta_tilde)
  aug_names <- setdiff(colnames(X), names(theta_tilde))
  
  #' Step 1: calculate W residuals
  mod <- lm(glue('cbind({glue_collapse(aug_names, sep = ",")}) ~ {glue_collapse(orig_names, sep = "+")}'), data = data.frame(X))
  gamma <- t(coef(mod))
  
  W <- resid(mod) 
  if(q == 1) {
    W <- matrix(W, dimnames = list(NULL, aug_names))
  }
  
  #' Step 2: fit the reduced glm model 
  theta.X <- drop(X[,orig_names, drop = FALSE] %*% theta_tilde)
  mod2 <- glm(Y ~ ., family = "binomial", data.frame(Y, theta.X, W))
  alpha <- coef(mod2)
  
  conv = (max(summary(mod2)$coefficients[,"Std. Error"]) < 1e3);
  #' Step 3: obtain final estimates
  intercept_hat <- alpha["(Intercept)"] - sum(alpha[aug_names] * gamma[,"(Intercept)"]) #intercept estimate
  
  beta_hat <- 
    c(alpha["theta.X"] * theta_tilde - 
        alpha[aug_names] %*% gamma[,orig_names],
      alpha[aug_names])
  names(beta_hat) <- c(orig_names, aug_names)
  
  return(list(intercept_hat = intercept_hat,
              beta_hat = beta_hat[colnames(X)],
              converge = conv, 
              message = NA,
              iter = NA,
              singular_hessian = NA,
              final_diff = NA, 
              objective_function = c("total" = NA)))
  
  
}
