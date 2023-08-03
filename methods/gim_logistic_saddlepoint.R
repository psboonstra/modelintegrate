#Pedro Orozco del Pino
#Using directly optim 

fxn_gim_logistic_saddlepoint <- function(Y, X_orig, X_aug,
                                         theta_tilde_with_intercept, beta_init_with_intercept, n_hist, 
                                         tol, max_rep, max_step_halfs = 10, max_lambda = Inf,
                                         Sigma0 = NULL) {
  
  ########################################### Define Sigma_0
  ## The GIM algorithm needs a matrix for the quadratic penalization
  # We can pre-calculate and fix Sig.hat.inv in some cases:
  if(!is.null(Sigma0) && is.matrix(Sigma0)) {
    #this approach is using the covariance of the historical data
    if(isSymmetric(Sigma0) && all(eigen(Sigma0, only.values = TRUE)$values > 0)) {
      Sig.hat.inv = solve(Sigma0) 
    }
    else {
      stop("May only provide symmetric, positive definite Sigma0")
    }
  } else if(!is.null(Sigma0) && Sigma0 == "likelihood") {
    #using the covariance of the historical model estimated with the current data 
    #remember X_orig_1 has a columns of ones
    # 12/21/22: Sig.hat.inv is missing an important multiplicative factor n
    Sig.hat.inv <- solve(vcov(glm(Y ~ -1 + X_orig_1, family = "binomial"))) / n
    #} else if(Sigma0=="true"&!is.null(beta.o)){
    # 12/21/22: I don't understand how beta.o comes into play here
    #using the true covariance of the estimate cov(hat.beta) = X' var(Y) X
    #Sigma0 <- t(X_orig_1)%*%diag(as.vector(prob_theta*(1-prob_theta)))%*%X_orig_1
    #V.inv <- solve(Sigma0)
  } else if(!is.null(Sigma0)) {
    stop("Wrong option of Sigma0: must be a fixed PD matrix, 'likelihood', or left as NULL for GIM estimation")
  }
  
  #The M.beta.theta function contains the pseudo likelihood to optimize  
  M.beta.theta <- function(beta_theta, n_hist = n_hist, return_value_only = TRUE) {
    #Remember that X_all_1 and X_orig_1 now contain a column of ones
    beta = beta_theta[1:p_all]
    theta = beta_theta[(p_all + 1):(p_all + p_orig)]
    linear_predictor = drop(X_all_1%*%beta);
    prob_beta = plogis(linear_predictor)
    prob_theta = plogis(drop(X_orig_1%*%theta))
    diff = drop(prob_beta - prob_theta)
    g.x.mu = X_orig_1 * diff
    
    # this is the likelihood of the current study ignoring external information
    loglikelihood = sum(Y * linear_predictor - log1plex(linear_predictor))
    
    #As proposed by Peisong run an inner loop to obtain a value for the lagrange multiplier
    #The Inn.loop function has to be run before any time we calculate the constraint of the CML approach
    foo <- Inner.loop(beta = beta, theta = theta)
    if(foo$singular) {
      stop("singular matrix in inner loop")
    }
    lambda = foo$lambda
    constraint = -sum(log1p(-g.x.mu%*%lambda))
    
    if(is.null(Sigma0)) {
      #this estimate is described in Zhang et al. 2020, is estimating the covariance by equations weighted with 
      #the estimate empirical distribution of the constrained optimization
      A =  A.hat(beta = beta,theta = theta,lambda = lambda) 
      B.inv = solve(B.hat(beta = beta,theta = theta,lambda = lambda))
      Sig.hat.inv = A %*% B.inv %*% A
    } 
    
    #######################################################
    #once the matrix of the quadratic penalization is set we then have a measured for the uncertainty of theta
    uncertainty.theta = -n_hist * t(theta_tilde_with_intercept - theta) %*% Sig.hat.inv %*% (theta_tilde_with_intercept - theta) / 2
    if(return_value_only) {
      return(loglikelihood + constraint + uncertainty.theta)
    } else {
      return(list(loglikelihood = loglikelihood, 
                  constraint = constraint, 
                  uncertainty.theta = uncertainty.theta,
                  lambda = lambda, 
                  Sig.hat.inv = Sig.hat.inv))
    }
  }
  #The inner loop is the function that obtains a value of lambda for a fixed beta and theta  
  Inner.loop <- function(beta, theta) {
    # See note below; I don't understand purpose of max_lambda
    # 12/22/22 but it seems to make a difference
    # max_lambda <- 2;
    
    prob_beta <- plogis(drop(X_all_1 %*% beta))
    prob_theta <- plogis(drop(X_orig_1 %*% theta))
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
      # Note also that you can still end up with lambda > max_lambda if 
      # the next while loop gives up when i == max_step_halfs
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
  Step1.outer.loop <- function(beta_theta, n_hist){
    m <- try(optim(beta_theta, M.beta.theta,
                   n_hist = n_hist,
                   control = list(fnscale = -1, maxit = 1),
                   method = "BFGS"), silent = TRUE)
    if("try-error" %in% class(m)) {
      m <- optim(beta_theta, M.beta.theta,
                 n_hist = n_hist,
                 control = list(fnscale = -1, maxit = 1),
                 method = "Nelder-Mead")
    }
    
    beta_theta_temp <- m$par
    delta <- beta_theta - beta_theta_temp
    return( delta )
  }
  #Step 2 checks is the step 1 step is valid if not then reduces the step length and tries again  
  Step2.outer.loop <- function(delta, beta_theta, n_hist){
    tau <- 1
    beta_theta_temp <- beta_theta - tau*delta
    M <- M.beta.theta(beta_theta, n_hist)
    M.temp <- M.beta.theta(beta_theta_temp, n_hist)
    k <- 1
    while(M.temp < M && k < max_step_halfs){
      tau <- tau / 2
      beta_theta_temp <- beta_theta - tau*delta
      M.temp <- M.beta.theta(beta_theta_temp, n_hist)
      k <- k + 1
    }
    return(beta_theta_temp)
  }
  #A hat is part one of the two parts of the estimation of the variance for the quadratic term, for details go to Zhang et al. 2020  
  A.hat<-function(beta, theta, lambda){
    prob_beta <- plogis(X_all_1%*%beta)
    prob_theta <- plogis(X_orig_1%*%theta)
    diff <- drop(prob_beta - prob_theta)
    g.x.mu <- X_orig_1 * diff
    
    p <- 1 / (n * (1 - g.x.mu%*%lambda))
    W <- as.vector(p*prob_theta)
    
    #The calculation notes have a minus but I got rid of it because it doesn't matter
    A <- t(X_orig_1)%*%diag(W)%*%X_orig_1
    
    return( as.matrix(A) )
  }
  #B hat is part two of the two parts of the estimation of the variance for the quadratic term, for details go to Zhang et al. 2020  
  B.hat<-function(beta, theta, lambda){
    prob_beta <- plogis(X_all_1%*%beta)
    prob_theta <- plogis(X_orig_1%*%theta)
    diff <- drop(prob_beta - prob_theta)
    g.x.mu <- X_orig_1 * diff
    p <- 1 / (n * (1 - g.x.mu%*%lambda))
    E.y2 <- prob_beta * (1 - prob_theta - prob_theta * (1 - prob_theta/prob_beta))
    W <- as.vector(p*E.y2)
    
    B <- t(X_orig_1)%*%diag(W)%*%X_orig_1
    
    return( as.matrix(B))
  }
  
  #Finally, we can actually begin the algorithm
  #extract the sample size of the current data
  n <- length(Y)
  
  #setting the union of both estimates as the initial point for the algorithm
  beta_theta_hat <- c(beta_init_with_intercept, theta_tilde_with_intercept)
  
  #Include a column of ones in the original covariates
  X_orig_1 <- cbind(1, X_orig)
  
  #Include a column of ones in the full data set of the current study
  X_all_1 <- cbind(X_orig_1, X_aug)
  
  # setting the dimensions of the full set of covariates and only the 
  # added covariates, in this sense both have the intercept so be aware of it
  p_orig <- ncol(X_orig_1)
  p_all <- ncol(X_all_1)
  length_theta_tilde_with_intercept <- length(theta_tilde_with_intercept)
  
  #initializing the errors
  max_change <- Inf

  #initializing the iteration counter
  iter <- 1
  singular <- F
  #n_hist_seq <- c(exp(seq(0, log(n_hist), length = min(max_rep, 25))), rep(n_hist, max_rep - min(max_rep, 25)) )
  
  while(tol < max_change  && iter <= max_rep){
    delta <- try(Step1.outer.loop(beta_theta_hat, n_hist = n_hist), silent = TRUE)
    if("try-error" %in% class(delta)) {
      if(str_detect(as.character(attributes(delta)$condition), "singular")) {
        singular <- TRUE
      } else {
        stop(delta[[1]])
      }
      max_change <- Inf
      break;
    }
    beta_theta_temp <- Step2.outer.loop(delta, beta_theta_hat, n_hist = n_hist)
    max_change <- max(abs(beta_theta_hat - beta_theta_temp))
    beta_theta_hat <- beta_theta_temp
    iter <- iter + 1
  }
  conv <- F
  if(max_change < tol && iter < max_rep && !singular){
    conv <- T
  }
  foo <- M.beta.theta(beta_theta_hat, n_hist = n_hist, return_value_only = FALSE)
  objective_function <- c("loglikelihood" = foo$loglikelihood,
                          "constraint" = foo$constraint,
                          "uncertainty.theta" = foo$uncertainty.theta, 
                          "total" = foo$loglikelihood + foo$constraint + foo$uncertainty.theta)
  return(list(intercept_hat = beta_theta_hat[1],
              beta_hat = beta_theta_hat[2:p_all],
              theta_hat_with_intercept = beta_theta_hat[(p_all+1):(p_all + p_orig)],
              beta_init_with_intercept = beta_init_with_intercept,
              theta_tilde_with_intercept = theta_tilde_with_intercept,
              lambda_hat = foo$lambda, 
              Sig.hat.inv = foo$Sig.hat.inv,
              objective_function = objective_function,
              iter = iter,
              converge = conv, 
              singular_hessian = singular,
              final_diff = max_change))
  
}

