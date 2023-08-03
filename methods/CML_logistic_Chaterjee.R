#Pedro Orozco del Pino 
#Code heavily based in code provided by Tian Gu 


###################### CML Algorithm
# is based on doing Newthon Raphson in the pseudo likelihood
# for more details follow S.2 of the Supplementary information of Chatterjee et al. 2016

fxn_CML_Logistic_Chat <- function(Y, X_orig, X_aug, 
                                  theta_hat, beta_hat_init, 
                                  tol, max_rep, step = 1) {
  
  n = length(Y)
  X_orig_1 = cbind(1, X_orig)
  X_all_i = cbind(X_orig_1, X_aug)
  p_orig <- ncol(X_orig_1)
  p_all <- ncol(X_all_i)
  beta_hat = beta_hat_init
  lambda = numeric(p_orig)
  estDiff = Inf
  counter = 0
  singular = F
  
  prob_theta_all = plogis(drop(X_orig_1 %*% theta_hat))
  
  while(estDiff > tol && counter < max_rep){
    v1 = 
      v2 = 
      h11 = 
      h12 = 
      h22 = 0
    
    for(i in 1:n) {
      
      X_orig_1_i = X_orig_1[i,]
      X_all_i_i = X_all_i[i,]
      Y_i = Y[i]
      
      prob_beta = plogis(sum(X_all_i_i * beta_hat))
      
      diff_prob = prob_beta - prob_theta_all[i]
      var_prob_beta = prob_beta * (1 - prob_beta)
      
      lam_T_X_orig_1_i = sum(X_orig_1_i * lambda)
      
      denom = as.numeric(1 - lam_T_X_orig_1_i * diff_prob)
      
      # score functions
      v1 = v1 + 
        # s_beta:
        (Y_i - prob_beta) * X_all_i_i +
        # tilde s_beta:
        lam_T_X_orig_1_i * var_prob_beta * X_all_i_i / denom
      
      v2 = v2 +
        # s_lambda
        diff_prob * X_orig_1_i / denom
      
      #hessian matrix
      h11 = h11 + 
        var_prob_beta * 
        (-1 + 
           (1 - 2 * prob_beta) * lam_T_X_orig_1_i / denom + 
           var_prob_beta * lam_T_X_orig_1_i^2 / denom^2) * 
        tcrossprod(X_all_i_i)
      
      # 12/19/22: line below should divide by denom^2 not denom (-phil)
      h12 = h12 + 
        #var_prob_beta * tcrossprod(X_orig_1_i, X_all_i_i) / denom
        var_prob_beta * tcrossprod(X_orig_1_i, X_all_i_i) / denom^2
      # Old way (equivalent)
      #(denom * var_prob_beta + var_prob_beta * diff_prob * lam_T_X_orig_1_i) * 
      #  tcrossprod(X_orig_1_i, X_all_i_i) / denom^2
      
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
    old_pars = c(beta_hat,lambda)
    new_pars = old_pars - step*drop(Hinv%*%V)
    estDiff = max(abs(new_pars - old_pars))
    beta_hat = new_pars[1:p_all]
    lambda = new_pars[(p_all+1):(p_all+p_orig)]
    counter = counter + 1
  }
  if(counter >= max_rep || singular ) { conv = F } else {conv = T}
  if(singular) {
    beta_hat = rep(NA, p_all)
  }
  list(intercept_hat = beta_hat[1],
       beta_hat = beta_hat[-1],
       beta_init = beta_hat_init[-1],
       theta_hat = theta_hat,
       iter = counter,
       converge = conv, 
       singular_hessian = singular,
       final_diff = estDiff)
  
}
