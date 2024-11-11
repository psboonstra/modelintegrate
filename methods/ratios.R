#' Ratios for Logistic Regression 
#' 
#' Implements the Ratios algorithm for logistic regression described in Taylor,
#' et al. (2022). Code is written by Pedro Orozco del Pino and Philip S.
#' Boonstra
#'
#'
#' @param Y vector of 0s and 1s
#' @param X design matrix with column names, with `nrow(X)` equal to `length(Y)`
#' @param theta_tilde named vector giving the constraint. Don't include the
#'   intercept. `names(theta_tilde)` should be a subset of  `colnames(X)`.
#'   
#' @return A named list. `intercept_hat` and `beta_hat` are estimated values of the
#' intercept and regression coefficients, respectively. 


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
