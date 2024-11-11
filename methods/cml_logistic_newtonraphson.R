#' CML for Logistic Regression using Newton-Raphson algorithms
#' 
#' Implements the CML algorithm for logistic regression using Newton-Raphson
#' algorithms as described in Section S.2 of Chatterjee et al. (2016). Code is
#' based on code provided by Tian Gu, then subsequently heavily modified by
#' Pedro Orozco del Pino and Philip S. Boonstra
#'
#' @param Y vector of 0s and 1s
#' @param X design matrix with column names, with `nrow(X)` equal to `length(Y)`
#' @param theta_tilde_with_intercept named vector giving the constraint,
#'   including the intercept. `names(theta_tilde_with_intercept)`
#'   should be a subset of  `colnames(X)`. The intercept should be
#'   named `(Intercept)`.
#' @param beta_init_with_intercept starting value to use for the algorithm.
#'   `names(beta_init_with_intercept)[-1]` should be equal to `colnames(X)`, with
#'   the first element of `names(beta_init_with_intercept)` being named `(Intercept)`.
#' @param tol Tolerance for convergence, e.g. $1e-5$
#' @param max_rep Maximum number of iterations of the outer loop before giving up
#' @param step Optional multiplier of the step size. Default is 1. Use a smaller
#'   value to take safer but less efficient steps
#'
#' @return A named list. `intercept_hat` and `beta_hat` are estimated values of the
#' intercept and regression coefficients, respectively. 


fxn_cml_logistic_newtonraphson <- function(Y, X,
                                           theta_tilde_with_intercept, beta_init_with_intercept, 
                                           tol, max_rep = 1000, step = 1) {
  
  stopifnot(all(colnames(X) == names(beta_init_with_intercept)[-1]))
  stopifnot(all(names(theta_tilde_with_intercept) %in% names(beta_init_with_intercept)))
  
  n <- length(Y)
  
  orig_var_names <- setdiff(names(theta_tilde_with_intercept), "(Intercept)")
  num_orig_plus1 <- length(theta_tilde_with_intercept)
  
  X_1 <- cbind(1, X)
  X_orig_1 <- cbind(1, X[, orig_var_names, drop = FALSE])
  
  num_all_plus1 <- ncol(X_1)
  beta_hat <- beta_init_with_intercept
  lambda_hat <- numeric(num_orig_plus1)
  max_change <- Inf
  iter <- 1
  singular <- F
  
  prob_theta <- plogis(drop(X_orig_1 %*% theta_tilde_with_intercept))
  
  pseudologlik <- matrix(0, nrow = max_rep, ncol = 2)
  message <- NA
  
  while(tol < max_change && iter <= max_rep) {
    
    v1 = 
      v2 = 
      h11 = 
      h12 = 
      h22 = 0
    
    for(i in 1:n) {
      
      X_orig_1_i = X_orig_1[i,]
      X_1_i = X_1[i,]
      Y_i = Y[i]
      
      linpred_i <- sum(X_1_i * beta_hat)
      prob_beta_i = plogis(linpred_i)
      pseudologlik[iter, 1] = 
        pseudologlik[iter, 1] + Y_i * linpred_i - log1plex(linpred_i)
      
      diff_prob = prob_beta_i - prob_theta[i]
      var_prob_beta_i = prob_beta_i * (1 - prob_beta_i)
      
      lam_T_X_orig_1_i = sum(X_orig_1_i * lambda_hat)
      
      denom = as.numeric(1 - lam_T_X_orig_1_i * diff_prob)
      pseudologlik[iter, 2] = 
        pseudologlik[iter, 2] - log(denom);
      
      # score functions
      v1 = v1 + 
        # s_beta:
        (Y_i - prob_beta_i) * X_1_i +
        # tilde s_beta:
        lam_T_X_orig_1_i * var_prob_beta_i * X_1_i / denom
      
      v2 = v2 +
        # s_lambda
        diff_prob * X_orig_1_i / denom
      
      #hessian matrix
      h11 = h11 + 
        var_prob_beta_i * 
        (-1 + 
           (1 - 2 * prob_beta_i) * lam_T_X_orig_1_i / denom + 
           var_prob_beta_i * lam_T_X_orig_1_i^2 / denom^2) * 
        tcrossprod(X_1_i)
      
      h12 = h12 + 
        var_prob_beta_i * tcrossprod(X_orig_1_i, X_1_i) / denom^2

      h22 = h22 + 
        diff_prob^2 * tcrossprod(X_orig_1_i) / denom^2 
      
    }
    #
    V = c(v1,v2)
    H = rbind(cbind(h11, t(h12)),cbind(h12, h22))
    Hinv <- try(solve(H), silent = T)
    if("try-error" %in% class(Hinv)) {
      singular = T
      message = Hinv[[1]]
      break;
    }
    old_pars = c(beta_hat, lambda_hat)
    new_pars = old_pars - step * drop(Hinv%*%V)
    max_change = max(abs(new_pars - old_pars))
    beta_hat = new_pars[1:num_all_plus1]
    lambda_hat = new_pars[(num_all_plus1 + 1):(num_all_plus1 + num_orig_plus1)]
    iter = iter + 1
  }
  conv <- F
  if(max_change < tol && iter < max_rep && !singular){
    conv <- T
  }
  if(singular) {
    beta_hat = rep(NA, num_all_plus1)
  }
  objective_function <- c("loglikelihood" = pseudologlik[iter - 1, 1],
                          "constraint" = pseudologlik[iter - 1, 2],
                          "total" = sum(pseudologlik[iter - 1, ]))
  list(intercept_hat = beta_hat[1],
       beta_hat = beta_hat[-1],
       beta_init_with_intercept = beta_init_with_intercept,
       theta_tilde_with_intercept = theta_tilde_with_intercept,
       lambda_hat = lambda_hat,
       objective_function = objective_function,
       iter = iter,
       converge = conv, 
       message = message,
       singular_hessian = singular,
       final_diff = max_change)
  
}
