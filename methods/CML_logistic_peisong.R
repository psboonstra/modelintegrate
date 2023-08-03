#Pedro Orozco del Pino
#Constrained Maximun Likelihood Method 
#implementation of Han and Lawless 2019, 
#EMPIRICAL LIKELIHOOD ESTIMATION USING AUXILIARY SUMMARY INFORMATION WITH DIFFERENT COVARIATE DISTRIBUTIONS

fxn_CML_Logistic_Han<-function(Y, X_orig, X_aug,
                               theta_hat, beta_hat_init,
                               tol, max_rep){
  
  #main function to optimize  
  M.beta<-function(beta,theta_hat,Y,X_orig_1,X_all_1){
    linear_predictor <- drop(X_all_1 %*% beta)
    prob_beta = plogis(linear_predictor)
    prob_theta = plogis(drop(X_orig_1 %*% theta_hat))
    diff = prob_beta - prob_theta
    g.x.mu = X_orig_1 * diff
    
    loglikelihood = sum(Y * linear_predictor - log1plex(linear_predictor))
    #
    foo = Inner.loop(beta,theta_hat,X_orig_1,X_all_1,n)
    if(foo$singular) {
      stop("singular matrix in inner loop")
    }
    lam = foo$lam
    constrain = -sum(log1p(-g.x.mu%*%lam))
    
    return( loglikelihood + constrain)
  }
  
  Inner.loop<-function(beta,theta_hat,X_orig_1,X_all_1,n){
    tol <- 1e-4
    steps <- 5
    l.max <- 2
    prob_beta <- plogis(drop(X_all_1 %*% beta))
    prob_theta <- plogis(drop(X_orig_1 %*% theta_hat))
    diff <- prob_beta - prob_theta
    g.x.mu <- X_orig_1 * diff
    
    lam <- numeric(length(theta_hat))
    k <- 1
    err <- Inf
    singular <- F
    
    while(err > tol && k <= steps){
      denom <- drop(1 - g.x.mu %*% lam)
      L <- -sum(log(denom))
      
      foo <- g.x.mu / denom
      J.l <- colSums(foo)#check recycling
      H.l <- crossprod(foo)
      
      if("try-error" %in% class(try(H.l.inv <- solve(H.l), silent = T))) {
        singular = T
        break;
      }
      delta = H.l.inv %*% J.l
      tau = 1
      lam.temp = lam - tau * delta
      cond1 = all( ( 1 - g.x.mu %*% lam.temp ) > 1/n ) 
      if(cond1){
        cond2 <- -sum(log1p(-g.x.mu %*% lam.temp)) < L
      } else{
        cond2 <- FALSE
      }
      cond3 = all(abs(lam.temp) < l.max )
      i= 1
      while( ( !cond2 | !cond3 )  && i < 50 ){
        tau = tau / 2
        lam.temp = lam - tau * delta
        cond1 = all( ( 1 - g.x.mu %*% lam.temp ) > 1/n ) 
        if(cond1){
          cond2 <- -sum(log1p(-g.x.mu %*% lam.temp)) < L
        }else{
          cond2 <- FALSE
        }
        cond3 = all(abs(lam.temp) < l.max )
        i <- i + 1 
      }
      err = sum(abs(lam - lam.temp))
      lam = lam.temp
      
      k = k + 1
      
    }
    
    return(list(lam = lam, singular = singular))    
  }
  
  Step1.outer.loop<-function(beta, theta_hat, Y, X_orig_1, X_all_1){
    m <- try(optim(beta, M.beta,
                   theta_hat = theta_hat, Y = Y, X_orig_1 = X_orig_1, X_all_1 = X_all_1,
                   control = list(fnscale = -1, maxit = 1),
                   method = "BFGS"), silent = TRUE)
    if("try-error" %in% class(m)) {
      m <- optim(beta, M.beta,
                 theta_hat = theta_hat, Y = Y, X_orig_1 = X_orig_1, X_all_1 = X_all_1,
                 control = list(fnscale = -1, maxit = 1),
                 method = "Nelder-Mead")
    }
    
    beta_temp <- m$par
    delta <- beta - beta_temp
    return( delta )
  }
  
  Step2.outer.loop <- function(delta,beta,theta_hat,Y,X_orig_1,X_all_1){
    attempts = 5
    tau = 1
    beta_temp = beta - tau*delta
    M = M.beta(beta,theta_hat,Y,X_orig_1,X_all_1)
    M_temp = M.beta(beta_temp,theta_hat,Y,X_orig_1,X_all_1)
    k <- 1
    while(M_temp < M && k < attempts){
      tau = tau / 2
      beta_temp = beta - tau * delta
      M_temp = M.beta(beta_temp,theta_hat,Y,X_orig_1,X_all_1)
      k = k + 1
    }
    return(beta_temp)
  }
  
  n <- length(Y)

  beta <- beta_hat_init
  X_orig_1 <- cbind(1, X_orig)
  X_all_1 <- cbind(1, X_orig, X_aug)
  
  err <- Inf
  iter <- 1
  singular <- F
  
  while(tol < err && iter <= max_rep){
    lam = Inner.loop(beta,theta_hat,X_orig_1,X_all_1,n)
    
    delta <- try(Step1.outer.loop(beta,theta_hat,Y,X_orig_1,X_all_1), silent = T)
    if("try-error" %in% class(delta)) {
      if(str_detect(as.character(attributes(delta)$condition), "singular")) {
        singular = T
      } else {
        stop(delta[[1]])
      }
      err <- Inf
      break;
    }
    
    beta_temp <- Step2.outer.loop(delta,beta,theta_hat,Y,X_orig_1,X_all_1)
    
    err <- sum( abs( beta - beta_temp ) )
    beta = beta_temp
    iter = iter + 1
  }
  conv = F
  if(err < tol && iter < max_rep && !singular){
    conv = T
  }
  return(list(intercept_hat = beta[1],
              beta_hat = beta[-1],
              beta_init = beta_hat_init[-1],
              theta_hat = theta_hat,
              iter = iter,
              converge = conv,
              singular_hessian = singular,
              final_diff = err))
  
}
