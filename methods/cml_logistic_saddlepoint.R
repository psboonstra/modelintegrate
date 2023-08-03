#Pedro Orozco del Pino
#Constrained Maximun Likelihood Method 
#implementation of Han and Lawless 2019, 
#EMPIRICAL LIKELIHOOD ESTIMATION USING AUXILIARY SUMMARY INFORMATION WITH DIFFERENT COVARIATE DISTRIBUTIONS

fxn_cml_logistic_saddlepoint <- function(Y, X_orig, X_aug,
                                         theta_tilde_with_intercept, beta_init_with_intercept,
                                         tol, max_rep, max_step_halfs = 10, max_lambda = 2) {
  
  #main function to optimize  
  M.beta <- function(beta, return_value_only = TRUE) {
    linear_predictor <- drop(X_all_1 %*% beta)
    prob_beta <- plogis(linear_predictor)
    diff <- prob_beta - prob_theta
    g.x.mu <- X_orig_1 * diff
    
    loglikelihood = sum(Y * linear_predictor - log1plex(linear_predictor))
    #
    foo = Inner.loop(beta)
    if(foo$singular) {
      stop("singular matrix in inner loop")
    }
    lambda = foo$lambda
    constraint = -sum(log1p(-g.x.mu%*%lambda))
    if(return_value_only) {
      return(loglikelihood + constraint)
    } else {
      return(list(loglikelihood = loglikelihood, 
                  constraint = constraint, 
                  lambda = lambda))
    }
    
  }
  
  Inner.loop <- function(beta) {
    # See note below; I don't understand purpose of max_lambda
    # max_lambda <- 2;
    prob_beta <- plogis(drop(X_all_1 %*% beta))
    diff <- prob_beta - prob_theta
    g.x.mu <- X_orig_1 * diff
    
    lambda <- numeric(length_theta_tilde_with_intercept)
    k <- 1
    max_change <- Inf
    singular <- F
    
    while(tol < max_change && k <= max(10, max_rep / 20)) {
      g.x.mu.lam <- drop(g.x.mu %*% lambda)
      denom <- 1 - g.x.mu.lam
      L <- -sum(log1p(-g.x.mu.lam))
      
      foo <- g.x.mu / denom
      J.l <- colSums(foo)#check recycling
      H.l <- crossprod(foo)
      
      if("try-error" %in% class(try(H.l.inv <- solve(H.l), silent = TRUE))) {
        singular <- T
        break;
      }
      delta <- H.l.inv %*% J.l
      tau <- 1
      lambda.temp <- lambda - tau * delta
      g.x.mu.lambda.temp <- g.x.mu %*% lambda.temp
      cond1 <- all( ( 1 - g.x.mu.lambda.temp ) > 1/n ) 
      if(cond1){
        cond2 <- -sum(log1p(-g.x.mu.lambda.temp)) < L
      } else{
        cond2 <- FALSE
      }
      # 12/21/22: I don't see condition 3 in the Inner Loop of Han and Lawless
      cond3 <- all(abs(lambda.temp) < max_lambda )
      i <- 1
      while( ( !cond2 | !cond3 ) && i < max_step_halfs) {
        tau <- tau / 2
        lambda.temp <- lambda - tau * delta
        g.x.mu.lambda.temp <- g.x.mu %*% lambda.temp
        cond1 <- all( ( 1 - g.x.mu.lambda.temp ) > 1/n ) 
        if(cond1){
          cond2 <- -sum(log1p(-g.x.mu.lambda.temp)) < L
        }else{
          cond2 <- FALSE
        }
        cond3 <- all(abs(lambda.temp) < max_lambda )
        i <- i + 1 
      }
      max_change <- sum(abs(lambda - lambda.temp))
      lambda <- lambda.temp
      
      k <- k + 1
    }
    
    return(list(lambda = lambda, singular = singular))    
  }
  
  #The step 1 of the outer loop find the step size   
  Step1.outer.loop <- function(beta) {
    m <- try(optim(beta, M.beta,
                   control = list(fnscale = -1, maxit = 1),
                   method = "BFGS"), silent = TRUE)
    if("try-error" %in% class(m)) {
      m <- optim(beta, M.beta,
                 control = list(fnscale = -1, maxit = 1),
                 method = "Nelder-Mead")
    }
    
    beta.temp <- m$par
    delta <- beta - beta.temp
    return( delta )
  }
  #Step 2 checks is the step 1 step is valid if not then reduces the step length and tries again  
  Step2.outer.loop <- function(delta, beta) {
    tau <- 1
    beta.temp <- beta - tau * delta
    M <- M.beta(beta)
    M.temp <- M.beta(beta.temp)
    k <- 1
    while(M.temp < M && k < max_step_halfs){
      tau <- tau / 2
      beta.temp <- beta - tau * delta
      M.temp <- M.beta(beta.temp)
      k <- k + 1
    }
    return(beta.temp)
  }
  
  n <- length(Y)
  
  X_orig_1 <- cbind(1, X_orig)
  X_all_1 <- cbind(1, X_orig, X_aug)
  beta_hat <- beta_init_with_intercept

  max_change <- Inf
  iter <- 1
  singular <- F
  prob_theta <- plogis(drop(X_orig_1 %*% theta_tilde_with_intercept))
  length_theta_tilde_with_intercept <- length(theta_tilde_with_intercept)
  
  while(tol < max_change && iter <= max_rep){
    
    delta <- try(Step1.outer.loop(beta_hat), silent = TRUE)
    if("try-error" %in% class(delta)) {
      if(str_detect(as.character(attributes(delta)$condition), "singular")) {
        singular <- T
      } else {
        stop(delta[[1]])
      }
      max_change <- Inf
      break;
    }
    
    beta_temp <- Step2.outer.loop(delta, beta_hat)
    
    max_change <- max(abs(beta_hat - beta_temp))
    beta_hat <- beta_temp
    iter <- iter + 1
  }
  conv <- F
  if(max_change < tol && iter < max_rep && !singular){
    conv <- T
  }
  foo <- M.beta(beta_hat, return_value_only = FALSE)
  objective_function <- c("loglikelihood" = foo$loglikelihood,
                          "constraint" = foo$constraint,
                          "total" = foo$loglikelihood + foo$constraint)
  return(list(intercept_hat = beta_hat[1],
              beta_hat = beta_hat[-1],
              beta_init_with_intercept = beta_init_with_intercept,
              theta_tilde_with_intercept = theta_tilde_with_intercept,
              lambda_hat = foo$lambda, 
              objective_function = objective_function,
              iter = iter,
              converge = conv,
              singular_hessian = singular,
              final_diff = max_change))
  
}
