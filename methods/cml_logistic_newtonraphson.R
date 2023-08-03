#Pedro Orozco del Pino 
#Code heavily based in code provided by Tian Gu 


###################### CML Algorithm
# is based on doing Newthon Raphson in the pseudo likelihood
# for more details follow S.2 of the Supplementary information of Chatterjee et al. 2016

fxn_cml_logistic_newtonraphson <- function(Y, X_orig, X_aug, 
                                           theta_tilde_with_intercept, beta_init_with_intercept, 
                                           tol, max_rep, step = 1) {
  
  n <- length(Y)
  X_orig_1 <- cbind(1, X_orig)
  X_all_1 <- cbind(X_orig_1, X_aug)
  p_orig <- ncol(X_orig_1)
  p_all <- ncol(X_all_1)
  beta_hat <- beta_init_with_intercept
  lambda_hat <- numeric(p_orig)
  max_change <- Inf
  iter <- 1
  singular <- F
  
  prob_theta_all <- plogis(drop(X_orig_1 %*% theta_tilde_with_intercept))
  
  pseudologlik2 <- matrix(0, nrow = max_rep, ncol = 2)
  
  while(tol < max_change && iter <= max_rep) {
    
    v1 = 
      v2 = 
      h11 = 
      h12 = 
      h22 = 0
    
    for(i in 1:n) {
      
      X_orig_1_i = X_orig_1[i,]
      X_all_1_i = X_all_1[i,]
      Y_i = Y[i]
      
      linpred <- sum(X_all_1_i * beta_hat)
      prob_beta = plogis(linpred)
      pseudologlik2[iter, 1] = 
        pseudologlik2[iter, 1] + Y_i * linpred - log1plex(linpred)
      
      diff_prob = prob_beta - prob_theta_all[i]
      var_prob_beta = prob_beta * (1 - prob_beta)
      
      lam_T_X_orig_1_i = sum(X_orig_1_i * lambda_hat)
      
      denom = as.numeric(1 - lam_T_X_orig_1_i * diff_prob)
      pseudologlik2[iter, 2] = 
        pseudologlik2[iter, 2] - log(denom);
      
      # score functions
      v1 = v1 + 
        # s_beta:
        (Y_i - prob_beta) * X_all_1_i +
        # tilde s_beta:
        lam_T_X_orig_1_i * var_prob_beta * X_all_1_i / denom
      
      v2 = v2 +
        # s_lambda
        diff_prob * X_orig_1_i / denom
      
      #hessian matrix
      h11 = h11 + 
        var_prob_beta * 
        (-1 + 
           (1 - 2 * prob_beta) * lam_T_X_orig_1_i / denom + 
           var_prob_beta * lam_T_X_orig_1_i^2 / denom^2) * 
        tcrossprod(X_all_1_i)
      
      # 12/19/22: line below should divide by denom^2 not denom (-phil)
      h12 = h12 + 
        #var_prob_beta * tcrossprod(X_orig_1_i, X_all_1_i) / denom
        var_prob_beta * tcrossprod(X_orig_1_i, X_all_1_i) / denom^2
      # Old way (equivalent)
      #(denom * var_prob_beta + var_prob_beta * diff_prob * lam_T_X_orig_1_i) * 
      #  tcrossprod(X_orig_1_i, X_all_1_i) / denom^2
      
      h22 = h22 + 
        diff_prob^2 * tcrossprod(X_orig_1_i) / denom^2 
      
    }
    #
    V = c(v1,v2)
    H = rbind(cbind(h11, t(h12)),cbind(h12, h22))
    if("try-error" %in% class(try(Hinv <- solve(H), silent = T))) {
      singular = T
      break;
    }
    old_pars = c(beta_hat,lambda_hat)
    new_pars = old_pars - step*drop(Hinv%*%V)
    max_change = max(abs(new_pars - old_pars))
    beta_hat = new_pars[1:p_all]
    lambda_hat = new_pars[(p_all+1):(p_all+p_orig)]
    iter = iter + 1
  }
  conv <- F
  if(max_change < tol && iter < max_rep && !singular){
    conv <- T
  }
  if(singular) {
    beta_hat = rep(NA, p_all)
  }
  objective_function <- c("loglikelihood" = pseudologlik2[iter-1, 1],
                          "constraint" = pseudologlik2[iter-1, 2],
                          "total" = sum(pseudologlik2[iter-1, ]))
  list(intercept_hat = beta_hat[1],
       beta_hat = beta_hat[-1],
       beta_init_with_intercept = beta_init_with_intercept,
       theta_tilde_with_intercept = theta_tilde_with_intercept,
       lambda_hat = lambda_hat,
       objective_function = objective_function,
       iter = iter,
       converge = conv, 
       singular_hessian = singular,
       final_diff = max_change)
  
}
