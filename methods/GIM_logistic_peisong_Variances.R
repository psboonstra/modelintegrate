#Pedro Orozco del Pino
#Using directly optim 

fxn_GIM_Logistic_Han <- function(Y, X_orig, X_aug,
                                 theta_hat, beta_hat_init, N, 
                                 tol, max_rep, Sigma0 = NULL, beta.o = NULL){
  #extract the sample size of the current data
  n <- length(Y)
  
  #setting the union of both estimates as the initial point for the algorithm
  beta.theta0 <- c(beta_hat_init, theta_hat)
  
  #I'm making a copy to have record of the initial point, this was thought to check is the method was moving at all
  beta.theta <- beta.theta0
  
  #Include a column of ones in the original covariates
  X_orig_1 <- cbind(1, X_orig)
  
  #Include a column of ones in the full data set of the current study
  X_all_1 <- cbind(X_orig_1, X_aug)
  
  # setting the dimensions of the full set of covariates and only the added covariates,
  # in this sense both have the intercept so be aware of it
  p_orig <- ncol(X_orig_1)
  p_all <- ncol(X_all_1)
  
  #initializing the errors
  err <- Inf
  err.M <- Inf
  
  #initializing the iteration counter
  iter = 1
  
  #The M.beta.theta function contains the pseudo likelihood to optimize  
  M.beta.theta <- function(beta.theta, theta_hat, Y, X_orig_1, X_all_1, N){
    #Remember that X_all_1 and X_orig_1 now contain a column of ones
    p_orig <- ncol(X_orig_1)
    p_all <- ncol(X_all_1)
    beta = beta.theta[1:p_all]
    theta = beta.theta[(p_all + 1):(p_all + p_orig)]
    linear_predictor = X_all_1%*%beta;
    Pb = plogis(linear_predictor)
    Pat = plogis(X_orig_1%*%theta)
    diff = drop(Pb - Pat)
    g.x.mu = X_orig_1 * diff
    
    # this is the likelihood of the current study ignoring external information
    loglikelihood = sum(Y*(linear_predictor) - log1plex(linear_predictor))
    
    # print(paste("Likelihood:",loglikelihood))
    #As proposed by Peisong run an inner loop to obtain a value for the lagrange multiplier
    foo <- Inner.loop(beta,theta,X_orig_1,X_all_1,n)
    if(foo$singular) {
      stop("singular matrix in inner loop")
    }
    lam = foo$lam
    #The Inn.loop function has to be run before any time we calculate the constrain of the CML approach
    constrain = -sum(log1p(-g.x.mu%*%lam))
    
    # print(paste("constrain:",constrain))
    
    ########################################### Define Sigma_0
    ## The GIM algorithm needs a matrix for the quadratic penalization, here are 4 options
    if(is.null(Sigma0)){
      #this estimate is described in Zhang et al. 2020, is estimating the covariance by equations weigthted with 
      #the estimate empirical distribution of the constrained optimization
      A.inv = solve( A.hat(X_orig_1,X_all_1,theta,beta,lam) )
      B = B.hat(X_orig_1,X_all_1,theta,beta,lam)
      Sig.hat = (A.inv%*%B%*%A.inv)
      
      V.inv = solve(Sig.hat)
    } else {
      if(!is.character(Sigma0)){
        #this approach is using the covariance of the historical data
        V.inv = solve(Sigma0)
      }else{
        if(Sigma0=="likelihood"){
          #using the covariance of the historical model estimated with the current data 
          Sigma0 <- vcov(glm(Y~X_orig_1-1,family="binomial"))#remember X_orig_1 has a columns of ones
          V.inv <- solve(Sigma0)
        }
        else{
          if(Sigma0=="true"&!is.null(beta.o)){
            #using the true covariance of the estimate cov(hat.beta) = X' var(Y) X
            Sigma0 <- t(X_orig_1)%*%diag(as.vector(Pat*(1-Pat)))%*%X_orig_1
            V.inv <- solve(Sigma0)
          }else{
            stop("Wrong option of Sigma0. Choose between true, likelihood or leave as NULL for GIM estimation")}
        }
      }
    }
    #######################################################
    #once the matrix of the quadratic penalization is set we then have a measured for the uncertainty of theta
    uncertainty.theta = -N*t(c(theta_hat-theta))%*%as.matrix(V.inv)%*%c(theta_hat-theta)/2
    # print(paste("uncer:",uncertainty.theta))
    value =  loglikelihood + constrain + uncertainty.theta
    
    return( value )
  }
  #The inner loop is the function that obtains a value of lambda for a fixed beta and theta  
  Inner.loop<-function(beta,theta,X_orig_1,X_all_1,n){
    tol = 1e-4
    steps = 10
    l.max <- 2
    
    Pb = plogis(X_all_1%*%beta)
    Pat = plogis(X_orig_1%*%theta)
    diff = drop(Pb - Pat)
    g.x.mu = X_orig_1 * diff
    
    lam = rep(0,length(theta))
    k = 1
    err = Inf
    singular = F
    while(err > tol && k <= steps) {
      denom = as.numeric( 1 - g.x.mu%*%lam )
      L = -sum( log(denom) )
      J.l = colSums( g.x.mu/denom )#check recycling
      H.l = crossprod(g.x.mu/denom)
      if("try-error" %in% class(try(H.l.inv <- solve(H.l), silent = T))) {
        singular = T
        break;
      }
      delta = H.l.inv %*% J.l
      
      tau = 1
      lam.temp = lam - tau*delta
      cond1 = all( ( 1 - g.x.mu%*%(lam.temp) ) > 1/n ) 
      if(cond1){
        cond2 <- -sum(log1p(-g.x.mu %*% lam.temp)) < L
      } else{
        cond2 <- FALSE
      }
      cond3 = all( abs(lam.temp) < l.max )
      i= 1
      while( ( !cond2 | !cond3 )  && i < 50 ){
        tau = tau / 2
        lam.temp = lam - tau*delta
        cond1 = all( ( 1 - g.x.mu%*%(lam.temp) ) > 1/n ) 
        if(cond1){
          cond2 <- -sum(log1p(-g.x.mu %*% lam.temp)) < L
        }else{
          cond2 <- FALSE
        }
        cond3 = all( abs(lam.temp) < l.max )
        i <- i + 1 
      }
      err = sum( abs( lam - lam.temp ) )
      lam = lam.temp
      k = k + 1
      
    }
    return(list(lam = lam, singular = singular))    
  }
  #The step 1 of the outer loop find the step size   
  Step1.outer.loop<-function(beta.theta,theta_hat,Y,X_orig_1,X_all_1,N){
    p = ncol(X_all_1)
    n = nrow(X_all_1)
    p_orig = ncol(X_orig_1)
    beta = beta.theta[1:p]
    theta = beta.theta[(p+1):(p+p_orig)]
    m <- try(optim(c(beta,theta), M.beta.theta,
                   theta_hat = theta_hat, Y=Y, X_orig_1=X_orig_1, X_all_1=X_all_1, N=N,
                   control = list(fnscale=-1,maxit=1),
                   method = "BFGS"), silent = TRUE)
    if("try-error" %in% class(m)) {
      m <- optim(c(beta,theta), M.beta.theta,
                 theta_hat = theta_hat, Y=Y, X_orig_1=X_orig_1, X_all_1=X_all_1, N=N,
                 control = list(fnscale=-1,maxit=1),
                 method = "Nelder-Mead")
    }
    
    beta.theta.temp <- m$par
    
    delta <- beta.theta - beta.theta.temp
    return( delta )
  }
  #Step 2 checks is the step 1 step is valid if not then reduces the step lenght and tries again  
  Step2.outer.loop<-function(delta,beta.theta,theta_hat,Y,X_orig_1,X_all_1,N){
    attempts = 5
    tau = 1
    beta.theta.temp = beta.theta -tau*delta
    M = M.beta.theta(beta.theta,theta_hat,Y,X_orig_1,X_all_1,N)
    M.temp = M.beta.theta(beta.theta.temp,theta_hat,Y,X_orig_1,X_all_1,N)
    k <- 1
    while(M.temp < M & k < attempts){
      tau = tau / 2
      beta.theta.temp = beta.theta -tau*delta
      M.temp = M.beta.theta(beta.theta.temp,theta_hat,Y,X_orig_1,X_all_1,N)
      k = k + 1
    }
    return(beta.theta.temp)
  }
  #A hat is part one of the two parts of the estimation of the variance for the quadratic term, for details go to Zhang et al. 2020  
  A.hat<-function(X_orig_1, X_all_1, theta, beta, lam){
    n = nrow(X_all_1)
    Pb = plogis(X_all_1%*%beta)
    Pat = plogis(X_orig_1%*%theta)
    diff = drop(Pb - Pat)
    g.x.mu = X_orig_1 * diff
    
    p = 1/(n*(1 - g.x.mu%*%lam))
    W = as.vector(p*Pat)
    
    #The calculation notes have a minus but I got rid of it because it doesn't matter
    A = t(X_orig_1)%*%diag(W)%*%X_orig_1
    
    return( as.matrix(A) )
  }
  #B hat is part two of the two parts of the estimation of the variance for the quadratic term, for details go to Zhang et al. 2020  
  B.hat<-function(X_orig_1,X_all_1,theta,beta,lam){
    n = nrow(X_all_1)
    Pb = plogis(X_all_1%*%beta)
    Pat = plogis(X_orig_1%*%theta)
    diff = drop(Pb - Pat)
    g.x.mu = X_orig_1 * diff
    p = 1/(n*(1 - g.x.mu%*%lam))
    E.y2 = Pb*(1-Pat-Pat*(1-Pat/Pb))
    W = as.vector(p*E.y2)
    
    B = t(X_orig_1)%*%diag(W)%*%X_orig_1
    
    return( as.matrix(B))
  }
  
  #Finally, we can actually begin the algorithm
  singular = F
  while(tol < err && iter <= max_rep ){
    beta = beta.theta[1:p_all]
    theta = beta.theta[(p_all+1):(p_all + p_orig)]
    delta <- try(Step1.outer.loop(c(beta,theta),theta_hat,Y,X_orig_1,X_all_1,N), silent = TRUE)
    if("try-error" %in% class(delta)) {
      if(str_detect(as.character(attributes(delta)$condition), "singular")) {
        singular = TRUE
      } else {
        stop(delta[[1]])
      }
      err <- Inf
      break;
    }
    beta.theta.temp <- Step2.outer.loop(delta,c(beta,theta),theta_hat,Y,X_orig_1,X_all_1,N)
    err <- sum( abs( beta.theta - beta.theta.temp ) )
    beta.theta = beta.theta.temp
    iter = iter + 1
  }
  conv = F
  if(err < tol && iter < max_rep && !singular){
    conv = T
  }
  return(list(intercept_hat = beta.theta[1],
              beta_hat = beta.theta[2:p_all],
              beta_init = beta_hat_init[-1],
              theta_hat = theta_hat[-1],
              theta_hat_updated = beta.theta[(p_all+1):(p_all + p_orig)],
              iter = iter,
              converge = conv, 
              singular_hessian = singular,
              final_diff = err))
  
}

