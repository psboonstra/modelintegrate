#Pedro Orozco del Pino
#Using directly optim 

fxn_gim_logistic_saddlepoint <- function(Y, X,
                                         theta_tilde_with_intercept, 
                                         beta_init_with_intercept, 
                                         theta_init_with_intercept = NULL, n_hist, 
                                         tol, max_rep = 1000, max_step_halfs = 25, max_abs_lambda = 2,
                                         Sigma_h = NULL, fix_beta = FALSE) {
  
  stopifnot(all(colnames(X)==names(beta_init_with_intercept)[-1]))
  stopifnot(all(names(theta_tilde_with_intercept) %in% names(beta_init_with_intercept)))
  if(is.null(theta_init_with_intercept)) {
    theta_init_with_intercept = theta_tilde_with_intercept
  } 
  stopifnot(all(names(theta_tilde_with_intercept) == names(theta_init_with_intercept)))
  
  ########################################### Define Sigma_0
  ## The GIM algorithm needs a matrix for the quadratic penalization
  # We can pre-calculate and fix chol_Sigma_h_inv in some cases:
  if(!is.null(Sigma_h) && is.matrix(Sigma_h)) {
    if(isSymmetric(Sigma_h, tol = sqrt(.Machine$double.eps), tol1 = NULL) && all(eigen(Sigma_h, only.values = TRUE)$values > 0)) {
      chol_Sigma_h_inv <- t(chol(solve(Sigma_h)))
    }
    else {
      stop("May only provide symmetric, positive definite Sigma_h")
    }
  } else if(!is.null(Sigma_h) && Sigma_h == "likelihood") {
    # Based upon the covariance using the current data
    orig_var_names <- setdiff(names(theta_tilde_with_intercept), "(Intercept)")
    X_orig <- X[, orig_var_names, drop = FALSE]
    chol_Sigma_h_inv <- t(chol(solve(vcov(glm(Y ~ X_orig, family = "binomial"))) / length(Y)))
    rm(orig_var_names, X_orig)
  } else if(!is.null(Sigma_h)) {
    stop("Wrong option of Sigma_h: must be a fixed PD matrix, 'likelihood', or left as NULL for GIM estimation")
  }
  
  #The M.beta.theta function contains the pseudo likelihood to optimize  
  M.beta.theta <- function(beta_theta, n_hist = n_hist, return_value_only = TRUE) {
    # X_1 and X_orig_1 contain a column of ones
    beta <- beta_theta[1:num_all_plus1]
    theta <- beta_theta[(num_all_plus1 + 1):(num_all_plus1 + num_orig_plus1)]
    linear_predictor <- drop(X_1%*%beta);
    prob_beta <- plogis(linear_predictor)
    prob_theta <- plogis(drop(X_orig_1%*%theta))
    diff <- drop(prob_beta - prob_theta)
    g.x.mu <- X_orig_1 * diff
    
    # this is the likelihood of the current study ignoring external information
    loglikelihood <- sum(Y * linear_predictor - log1plex(linear_predictor))
    
    # As proposed in Han and Lawless (2019), run an inner loop to obtain a value
    # for the Lagrange multiplier.
    foo <- Inner.loop(g.x.mu)
    if(foo$singular) {
      stop("singular matrix in inner loop")
    }
    lambda <- foo$lambda
    g.x.mu.lam <- drop(g.x.mu %*% lambda)
    constraint <- -sum(log1p(-g.x.mu.lam))
    
    if(is.null(Sigma_h)) {
      # this estimate is described in Zhang et al. 2020, is estimating the
      # covariance by equations weighted with the estimate empirical
      # distribution of the constrained optimization
      lagrange <- 1 / (n * (1 - g.x.mu.lam))
      chol_Sigma_h_inv <- make_half_sandwich(prob_beta, prob_theta, lagrange)
    } 
    
    # once the matrix of the quadratic penalization is set we then have a
    # measured for the uncertainty of theta
    quadratic.term <- (-n_hist) * tcrossprod(crossprod(theta_tilde_with_intercept - theta, chol_Sigma_h_inv)) / 2
    if(return_value_only) {
      return(loglikelihood + constraint + quadratic.term)
    } else {
      return(list(loglikelihood = loglikelihood, 
                  constraint = constraint, 
                  quadratic.term = quadratic.term,
                  lambda = lambda, 
                  chol_Sigma_h_inv = chol_Sigma_h_inv))
    }
  }
  # The inner loop is the function that obtains a value of lambda for a fixed
  # beta and theta. See Han and Lawless (2019)
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
  # Step 1 of the outer loop finds the step size. We reverse engineer the step
  # by calculating one step of a numerical solver and treating the difference
  # between the returned parameters and the initial parameters as the value
  # of delta
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
    
    delta <- (beta_theta - m$par);
    return( delta )
  }
  # Step 2 checks if the step 1 step is valid if not then reduces the step
  # length and tries again
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
  # This returns sandwich estimate of sigma^{-1}
  # For scrumptious details on the sandwich see Zhang et al. 2020
  make_half_sandwich <- function(prob_beta, prob_theta, lagrange){
    W_A <- lagrange * prob_theta * (1 - prob_theta)
    W_B <- lagrange * (prob_theta^2 + prob_beta - 2 * prob_theta * prob_beta)
    A <- crossprod(sqrt(W_A) * X_orig_1)
    B <- crossprod(sqrt(W_B) * X_orig_1)
    # This is A %*% B^{-1/2}
    A %*% backsolve(chol(B),identity_num_orig_plus1)
  }
  
  # Finally, we can actually begin the algorithm
  # extract the sample size of the current data
  n <- length(Y)
  
  # concatenate initial values 
  beta_theta_hat <- c(beta_init_with_intercept, theta_init_with_intercept)
  
  # Include a column of ones in the full data set of the current study
  # and the original covariates only
  X_1 <- cbind(1, X)
  orig_var_names <- setdiff(names(theta_tilde_with_intercept), "(Intercept)")
  X_orig_1 <- cbind(1, X[, orig_var_names, drop = FALSE])
  
  # setting the dimensions of the full set of covariates and only the 
  # added covariates, in this sense both have the intercept so be aware of it
  num_orig_plus1 <- ncol(X_orig_1)
  num_all_plus1 <- ncol(X_1)
  identity_num_orig_plus1 <- diag(1, num_orig_plus1)
  
  #initializing the errors
  max_change <- Inf
  
  #initializing the iteration counter
  iter <- 1
  singular <- F
  message <- NA
  
  while(tol < max_change  && iter <= max_rep){
    delta <- try(Step1.outer.loop(beta_theta_hat, n_hist = n_hist), silent = TRUE)
    if("try-error" %in% class(delta)) {
      if(str_detect(as.character(attributes(delta)$condition), "singular")) {
        singular <- TRUE
        message <- delta[[1]]
      } else {
        stop(delta[[1]])
      }
      max_change <- Inf
      break;
    }
    beta_theta_temp <- Step2.outer.loop(delta, beta_theta_hat, n_hist = n_hist)
    max_change <- max(abs(beta_theta_hat - beta_theta_temp))
    beta_theta_hat <- beta_theta_temp
    if(fix_beta) {
      beta_theta_hat[1:num_all_plus1] <- beta_init_with_intercept
    }
    iter <- iter + 1
  }
  conv <- F
  if(max_change < tol && iter < max_rep && !singular){
    conv <- T
  }
  if(conv) {
    foo <- M.beta.theta(beta_theta_hat, n_hist = n_hist, return_value_only = FALSE)
    objective_function <- c("loglikelihood" = foo$loglikelihood,
                            "constraint" = foo$constraint,
                            "quadratic.term" = foo$quadratic.term, 
                            "total" = foo$loglikelihood + foo$constraint + foo$quadratic.term)
    return(list(intercept_hat = beta_theta_hat[1],
                beta_hat = beta_theta_hat[2:num_all_plus1],
                theta_hat_with_intercept = beta_theta_hat[(num_all_plus1+1):(num_all_plus1 + num_orig_plus1)],
                beta_init_with_intercept = beta_init_with_intercept,
                theta_tilde_with_intercept = theta_tilde_with_intercept,
                lambda_hat = foo$lambda[,1], 
                Sigma_h_inv = tcrossprod(foo$chol_Sigma_h_inv),
                objective_function = objective_function,
                converge = conv, 
                message = message,
                iter = iter,
                singular_hessian = singular,
                final_diff = max_change))
  } else {
    return(list(intercept_hat = beta_theta_hat[1],
                beta_hat = beta_theta_hat[2:num_all_plus1],
                theta_hat_with_intercept = beta_theta_hat[(num_all_plus1+1):(num_all_plus1 + num_orig_plus1)],
                beta_init_with_intercept = beta_init_with_intercept,
                theta_tilde_with_intercept = theta_tilde_with_intercept,
                lambda_hat = NA, 
                Sigma_h_inv = NA,
                objective_function = NA,
                converge = conv, 
                message = message,
                iter = iter,
                singular_hessian = singular,
                final_diff = max_change))
  }
  
}


