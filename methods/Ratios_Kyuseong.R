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
#'theta_hat_no_intercept effect size estimates from the historical study for the original covariates (drop the intercept)
#' alpha_1 = 1 does not estimates alpha_1 but assumes it is equal to 1. In this case we include an offset in the estimation of alpha
#' Center X_orig


ratios <- function(Y, X_orig, X_aug, theta_hat_no_intercept, alpha_1=1) {
  
  X_orig_centered <- scale(X_orig, scale = FALSE)
  
  q <- dim(X_aug)[2]
  
  n_curr <- dim(X_aug)[1]
  p <- dim(X_orig)[2]
  
  #' calculate residuals
  W <- matrix(NA, ncol = q, nrow = n_curr)
  gamma <- matrix(NA, ncol = p + 1, nrow = q)
  for( i in 1:q){
    mod <- lm(X_aug[,i] ~ X_orig_centered)
    W[,i] <- mod$residuals
    gamma[i,] <- mod$coefficients
  }
  
  #' fit the glm model 

  theta.X <- X_orig_centered%*%theta_hat_no_intercept
  
  if(alpha_1 == 1) {
    fit <- glm(Y ~ W, offset = theta.X, family = "binomial");
  } else{
    fit <- glm(Y ~ theta.X + W, family = "binomial")
  }

  alpha <- fit$coefficients
  
  conv = (max(summary(fit)$coefficients[,2]) < 1e3);
  
  #' obtain final estimates
  
  if(alpha_1==1){
    beta0 <- alpha[1] - t(alpha[-c(1)])%*%gamma[,1] #intercept estimate
    beta <- c(NULL)
    for(j in 1:p){
      beta[j] <- theta_hat_no_intercept[j] - t(alpha[-c(1)])%*%gamma[,(j+1)]  #original covariates estimates
    }
    beta[(p+1):(p+q)] <- alpha[-c(1)] # added covariates estimates
  }else{
    beta0 <- alpha[1] - t(alpha[-c(1,2)])%*%gamma[,1] #intercept estimate
    beta <- c(NULL)
    for(j in 1:p){
      beta[j] <- alpha[2]*theta_hat_no_intercept[j] - t(alpha[-c(1,2)])%*%gamma[,(j+1)]  #original covariates estimates
    }
    beta[(p+1):(p+q)] <- alpha[-c(1,2)] # added covariates estimates
  }
 
  return(list(intercept_hat = drop(beta0),
              beta_hat = beta,
              converge = conv, 
              iter = NA,
              singular_hessian = NA,
              final_diff = NA))

  
}
