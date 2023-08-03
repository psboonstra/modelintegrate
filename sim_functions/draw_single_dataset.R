
# input: 
# 1. true_betas ( vector ) concatenated original + added covariates effect sizes
# 2. true_betas_hist ( vector ) concatenated original+added covariates effect 
# size for historical model (can be equal to true betas or not but must be
# same length)
# 3. num_orig, num_aug, num_all_coef (all scalars) these are integers
# that should satisfy num_orig + num_aug == num_all_coefs == length(true_betas)
# 4. n_curr (scalar) sample size of current study
# 5. n_hist (scalar) sample size of historical study
# 6. true_chol_sigma_x (matrix) cholesky decomposition of covariance matrix 
# of covariates from current/internal population
# 7. true_chol_sigma_x_hist (scalar) cholesky decomposition of covariance matrix 
# of covariates from historical/external population (can be equal to 
# true_chol_sigma_x or not)

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


Xc <-
  mvrnorm(n_curr, numeric(num_all_coef), diag(1, num_all_coef)) %*% true_chol_sigma_x

#probability of being a case conditional on covariates
Pc <- drop(plogis(true_beta0 + Xc%*%true_betas))

#Simulate the binary response using the vector of conditional probabilities
Yc <- rbinom(n_curr, 1, Pc)

#data of current study
Xc_orig <- Xc[, 1:num_orig, drop = F]#original covariates
Xc_aug <- Xc[, num_orig + (1:num_aug), drop = F]#added covariates

# This is the set of historical covariates for all 
Xh <-
  mvrnorm(n_hist, numeric(num_all_coef), diag(1, num_all_coef)) %*% true_chol_sigma_x_hist

# probability of being a case conditional on covariates
Ph <- drop(plogis(true_beta0 + Xh %*% true_betas_hist))

# binary response for the historical study
Yh <- rbinom(n_hist, 1, Ph)

# data of historical study
Xh_orig <- Xh[seq_len(n_hist), 1:num_orig, drop = F]#original covariates
Xh_aug <- Xh[seq_len(n_hist), num_orig + (1:num_aug), drop = F]#added covariates

mod_hist <- logistf(Yh ~ Xh_orig, family="binomial", plconf = 1)

# saving the estimated effect sizes, including the intercept
theta_tilde_with_intercept <- mod_hist$coefficients
theta_tilde <- theta_tilde_with_intercept[-1]

# saving the covariance matrix of the estimated effect sizes of the historical study
var_theta_tilde_with_intercept <- vcov(mod_hist)
var_theta_tilde <- var_theta_tilde_with_intercept[-1, -1, drop = FALSE]

# starting values for beta
beta_init_with_intercept <- coef(logistf(Yc ~ Xc_orig + Xc_aug,family="binomial", plconf = 1))


Xv <-
  mvrnorm(2.5e3, numeric(num_all_coef), diag(1, num_all_coef)) %*% true_chol_sigma_x

#probability of being a case conditional on covariates
Pv <- drop(plogis(true_beta0 + Xv%*%true_betas))

Yv <- rbinom(length(Pv), 1, Pv)
