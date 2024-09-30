#Pedro Orozco del Pino
#Constrained Maximun Likelihood Method 
#implementation of Han and Lawless 2019, 
#EMPIRICAL LIKELIHOOD ESTIMATION USING AUXILIARY SUMMARY INFORMATION WITH DIFFERENT COVARIATE DISTRIBUTIONS


fxn_cml_logistic_saddlepoint <- function(Y, X,
                                         theta_tilde_with_intercept, beta_init_with_intercept,
                                         tol, max_rep = 1000, max_step_halfs = 25, max_abs_lambda = 2) {
  
  stopifnot(all(colnames(X) == names(beta_init_with_intercept)[-1]))
  stopifnot(all(names(theta_tilde_with_intercept) %in% names(beta_init_with_intercept)))
  
  #main function to optimize  
  M.beta <- function(beta, return_value_only = TRUE) {
    linear_predictor <- drop(X_1 %*% beta)
    prob_beta <- plogis(linear_predictor)
    diff <- prob_beta - prob_theta
    g.x.mu <- X_orig_1 * diff
    
    # this is the likelihood of the current study ignoring external information
    loglikelihood = sum(Y * linear_predictor - log1plex(linear_predictor))
    
    # As proposed in Han and Lawless (2019), run an inner loop to obtain a value
    # for the Lagrange multiplier.
    foo = Inner.loop(g.x.mu)
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
  
  Inner.loop <- function(g.x.mu) {
    
    lambda <- numeric(num_orig_plus1)
    k <- 1
    max_change <- Inf
    singular <- F
    
    while(tol < max_change && k <= max(10, max_rep / 20)) {
      g.x.mu.lam <- drop(g.x.mu %*% lambda)
      denom <- 1 - g.x.mu.lam
      L <- -sum(log1p(-g.x.mu.lam))
      
      foo <- g.x.mu / denom
      J.l <- colSums(foo)
      H.l <- crossprod(foo)
      
      if("try-error" %in% class(try(H.l.inv <- chol2inv(chol(H.l)), silent = TRUE))) {
        singular <- T
        break;
      }
      delta <- H.l.inv %*% J.l
      tau <- 1
      lambda.temp <- lambda - tau * delta
      cond3 <- all(abs(lambda.temp) < max_abs_lambda )
      if(cond3) {
        g.x.mu.lambda.temp <- g.x.mu %*% lambda.temp
        cond1 <- all( ( 1 - g.x.mu.lambda.temp ) > 1/n ) 
        if(cond1){
          cond2 <- -sum(log1p(-g.x.mu.lambda.temp)) < L
        } else{
          cond2 <- FALSE
        } 
      } else {
        cond2 <- FALSE
      }
      i <- 1
      while( ( !cond2 | !cond3 ) && i < max_step_halfs) {
        tau <- tau / 2
        lambda.temp <- lambda - tau * delta
        i <- i + 1 
        cond3 <- all(abs(lambda.temp) < max_abs_lambda )
        if(!cond3) {next;}
        g.x.mu.lambda.temp <- g.x.mu %*% lambda.temp
        cond1 <- all( ( 1 - g.x.mu.lambda.temp ) > 1/n ) 
        if(cond1) {
          cond2 <- -sum(log1p(-g.x.mu.lambda.temp)) < L
        } else{
          cond2 <- FALSE
        }
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
  
  orig_var_names <- setdiff(names(theta_tilde_with_intercept), "(Intercept)")
  num_orig_plus1 <- length(theta_tilde_with_intercept)
  
  X_1 <- cbind(1, X)
  X_orig_1 <- cbind(1, X[, orig_var_names, drop = FALSE])
  
  beta_hat <- beta_init_with_intercept
  
  max_change <- Inf
  iter <- 1
  singular <- F
  prob_theta <- plogis(drop(X_orig_1 %*% theta_tilde_with_intercept))
  message <- NA
  
  while(tol < max_change && iter <= max_rep){
    
    delta <- try(Step1.outer.loop(beta_hat), silent = TRUE)
    if("try-error" %in% class(delta)) {
      if(str_detect(as.character(attributes(delta)$condition), "singular")) {
        singular <- T
        message <- delta[[1]]
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
              lambda_hat = foo$lambda[,1], 
              converge = conv,
              message = message, 
              iter = iter,
              singular_hessian = singular,
              final_diff = max_change,
              objective_function = objective_function))
  
}



