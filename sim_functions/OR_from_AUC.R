
auc_to_betas <- function(target_auc = 0.7, relative_scales, cor, n_samps = 3e6, seed = sample(.Machine$integer.max, 1)) {
  require(mnormt);
  #starting value
  set.seed(seed);
  X <- MASS::mvrnorm(n_samps, numeric(length(relative_scales)), cor + diag(1 - cor, length(relative_scales))) 
  log_beta_raw = -2;
  lin_pred <- drop(X %*% (exp(log_beta_raw) * relative_scales))
  P <- drop(plogis(lin_pred))
  Y <- rbinom(n_samps, 1, P)
  sample_events = sample(which(Y==1), round(n_samps / 2), replace = T)
  sample_non_events = sample(which(Y==0), round(n_samps / 2), replace = T)
  
  current_auc <- 
    mean(lin_pred[sample_events] > lin_pred[sample_non_events] + 1e-8) + 
    0.5 * mean(abs(lin_pred[sample_events] - lin_pred[sample_non_events]) < 1e-8) 
  
  i = 1;
  while(i < 100 && abs(log(target_auc / current_auc)) > 1e-4) {
    log_beta_raw <- log_beta_raw + log(target_auc / current_auc)
    
    lin_pred <- drop(X %*% (exp(log_beta_raw) * relative_scales))
    P <- drop(plogis(lin_pred))
    Y <- rbinom(n_samps, 1, P)
    sample_events = sample(which(Y==1), round(n_samps / 2), replace = T)
    sample_non_events = sample(which(Y==0), round(n_samps / 2), replace = T)
    
    current_auc <- 
      mean(lin_pred[sample_events] > lin_pred[sample_non_events] + 1e-8) + 
      0.5 * mean(abs(lin_pred[sample_events] - lin_pred[sample_non_events]) < 1e-8) 
    i = i + 1;
  }
  list(betas = round(exp(log_beta_raw) * relative_scales, 3),
       auc_est = current_auc)
}

betas_to_auc <- function(betas, cor, n_samps = 3e6, seed = sample(.Machine$integer.max, 1)) {
  require(mnormt);
  
  set.seed(seed);
  X <- MASS::mvrnorm(n_samps, numeric(length(betas)), cor + diag(1 - cor, length(betas))) 
  
  lin_pred <- drop(X %*% betas)
  P <- drop(plogis(lin_pred))
  Y <- rbinom(n_samps, 1, P)
  sample_events = sample(which(Y==1), round(n_samps / 2), replace = T)
  sample_non_events = sample(which(Y==0), round(n_samps / 2), replace = T)
  
  current_auc <- 
    mean(lin_pred[sample_events] > lin_pred[sample_non_events]) + 
    0.5 * mean(abs(lin_pred[sample_events] - lin_pred[sample_non_events]) < 1e-8) 
  
  current_auc
}


# Includes the intercept
betas_to_auc2 <- function(betas, cor, intercept = 0, n_samps = 3e6, seed = sample(.Machine$integer.max, 1)) {
  require(mnormt);
  
  set.seed(seed);
  X <- MASS::mvrnorm(n_samps, numeric(length(betas)), cor + diag(1 - cor, length(betas))) 
  
  lin_pred <- intercept + drop(X %*% betas)
  P <- drop(plogis(lin_pred))
  Y <- rbinom(n_samps, 1, P)
  sample_events = sample(which(Y==1), round(n_samps / 2), replace = T)
  sample_non_events = sample(which(Y==0), round(n_samps / 2), replace = T)
  
  current_auc <- 
    mean(lin_pred[sample_events] > lin_pred[sample_non_events] + 1e-8) + 
    0.5 * mean(abs(lin_pred[sample_events] - lin_pred[sample_non_events]) < 1e-8) 
  current_auc
}


betas_to_prevalence <- function(betas, cor, intercept = 0, n_samps = 3e6, seed = sample(.Machine$integer.max, 1)) {
  require(mnormt);
  
  set.seed(seed);
  X <- MASS::mvrnorm(n_samps, numeric(length(betas)), cor + diag(1 - cor, length(betas))) 
  
  lin_pred <- intercept + drop(X %*% betas)
  P <- drop(plogis(lin_pred))
  Y <- rbinom(n_samps, 1, P)
  mean(Y);
}



