
# input: 
# 1. true_betas ( vector ) concatenated original + added covariates effect sizes
# 2. true_betas_hist ( vector ) concatenated original+added covariates effect 
# size for historical model (can be equal to true betas or not but must be
# same length)
# 3. num_orig, num_aug, num_all (all scalars) these are integers
# that should satisfy num_orig + num_aug == num_all == length(true_betas)
# 4. n_curr (scalar) sample size of current study
# 5. n_hist (scalar) sample size of historical study
# 6. true_chol_sigma_x (matrix) cholesky decomposition of covariance matrix 
# of covariates from current/internal population
# 7. true_chol_sigma_x_hist (matrix) cholesky decomposition of covariance matrix 
# of covariates from historical/external population (can be equal to 
# true_chol_sigma_x or not)
# 8. all_var_names: character vector as long as num_orig+ num_aug
# to name the covariates. Should be distinct. 
# 9. which_binary: logical vector with length num_all indicating which
# covariates should be transformed to 0/2 values using the function 2*I[x>0]. The
# scale of 2 is used to make the resulting variance be approximately equal to 1

#output:
# 1. Xc_orig (matrix) simulated current original covariates
# 2. Xc_aug (matrix) simulated current added covariates
# 3. Yc (vector) simulated current outcome
# 4. Xh_orig (matrix) simulated historical original covariates
# 5. Xh_aug (matrix) simulated historical added covariates
# 6. Yh (vector) simulated historical outcome
# 7. theta_tilde (vector) historical effect size estimates with the historical model
# 8. var_theta_tilde (matrix) covariance matrix of theta_tilde. Note this is 
# assumed to be the limiting covariance matrix *divided by* n_hist. 


Xc <- (mvrnorm(n_curr, numeric(num_all), diag(1, num_all)) %*% true_chol_sigma_x) %>%
  `colnames<-`(all_var_names)
Xc[,which_binary] <- 2 * (Xc[,which_binary] > 0) 

#probability of being a case conditional on covariates
Pc <- drop(plogis(true_beta0 + Xc%*%true_betas))

#Simulate the binary response using the vector of conditional probabilities
Yc <- rbinom(n_curr, 1, Pc)

#data of current study
Xc_orig <- Xc[, 1:num_orig, drop = F]
Xc_aug <- Xc[, num_orig + (1:num_aug), drop = F]

# This is the set of historical covariates for all 
Xh <- (mvrnorm(n_hist, numeric(num_all), diag(1, num_all)) %*% true_chol_sigma_x_hist) %>%
  `colnames<-`(all_var_names)
Xh[,which_binary] <- 2 * (Xh[,which_binary] > 0) 

# probability of being a case conditional on covariates
Ph <- drop(plogis(true_beta0 + Xh %*% true_betas_hist))

# binary response for the historical study
Yh <- rbinom(n_hist, 1, Ph)

# data of historical study
Xh_orig <- Xh[, 1:num_orig, drop = F]
Xh_aug <- Xh[, num_orig + (1:num_aug), drop = F]

mod_hist <- logistf(Y ~ ., family="binomial", plconf = 1, data = data.frame(Y = Yh, Xh_orig))

# saving the estimated effect sizes, including the intercept
theta_tilde_with_intercept <- mod_hist$coefficients
theta_tilde <- theta_tilde_with_intercept[-1]

# saving the covariance matrix of the estimated effect sizes of the historical study
var_theta_tilde_with_intercept <- vcov(mod_hist)
var_theta_tilde <- var_theta_tilde_with_intercept[-1, -1, drop = FALSE]

# starting values for beta
beta_init_with_intercept <- coef(logistf(Y ~ ., family="binomial", plconf = 1, data = data.frame(Y = Yc, Xc)))

Xv <- (mvrnorm(2.5e3, numeric(num_all), diag(1, num_all)) %*% true_chol_sigma_x) %>%
  `colnames<-`(all_var_names)
Xv[,which_binary] <- 2 * (Xv[,which_binary] > 0) 

#probability of being a case conditional on covariates
Pv <- drop(plogis(true_beta0 + Xv%*%true_betas))

Yv <- rbinom(length(Pv), 1, Pv)
