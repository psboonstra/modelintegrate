###########Pedro Orozco del Pino
### Simulations for adaptive Bayes

fxn_Adapt_Bayes<-function(Y, X_orig, X_aug, 
                          theta_hat_no_intercept, #theta without the intercept
                          var_theta_hat_no_intercept, #var(theta) without the intercept
                          beta_scale,
                          phi_mean = 1,
                          phi_sd = 0,
                          mc_warmup = 1e3,
                          mc_iter_after_warmup = 1e3,
                          mc_chains = 2,
                          mc_thin = 1,
                          mc_stepsize = 0.1,
                          mc_adapt_delta = 0.999,
                          mc_max_treedepth = 15,
                          true_cor = NULL,
                          family = "binomial", 
                          stan_seed = sample(.Machine$integer.max, 1)){
  p = dim(X_orig)[2]
  q = dim(X_aug)[2]
  eigendecomp_hist_var = eigen(var_theta_hat_no_intercept);
  colnames(X_orig) = as.character(glue("p{1:ncol(X_orig)}"))
  colnames(X_aug) = as.character(glue("q{1:ncol(X_aug)}"))
  if(is.null(true_cor)){
    aug_projection <- 
      create_projection(x_curr_orig = X_orig,
                        x_curr_aug = X_aug,
                        imputes_list = list(c(1, 500)),
                        eigenvec_hist_var = t(eigendecomp_hist_var$vectors))[[1]]
  } else {
    E_xa.xo = true_cor / (1 - true_cor + true_cor * p)
    aug_projection <- matrix(E_xa.xo,ncol=q,nrow=p)
  }  
  scale_to_variance225 = diag(var_theta_hat_no_intercept) / 225; #What is that 225?
  mod <- glm_sab(y = Y,
                 x_standardized = cbind(X_orig, X_aug),
                 alpha_prior_mean = theta_hat_no_intercept, 
                 alpha_prior_cov = var_theta_hat_no_intercept,
                 aug_projection = aug_projection,
                 phi_mean = phi_mean,
                 phi_sd = phi_sd,
                 beta_orig_scale = beta_scale, 
                 beta_aug_scale = beta_scale,
                 local_dof = 1, 
                 global_dof = 1, 
                 slab_precision = (1/15)^2,
                 only_prior = F, 
                 mc_warmup = mc_warmup,
                 mc_iter_after_warmup = mc_iter_after_warmup,
                 mc_chains = mc_chains,
                 mc_thin = mc_thin,
                 mc_stepsize = mc_stepsize,
                 mc_adapt_delta = mc_adapt_delta,
                 mc_max_treedepth = mc_max_treedepth,
                 phi_dist = "trunc_norm",
                 family = family,
                 eigendecomp_hist_var = eigendecomp_hist_var,
                 scale_to_variance225 = scale_to_variance225, 
                 seed = stan_seed)
  
  beta_hat <- colMeans(mod$curr_beta)
  intercept_hat <- mean(mod$curr_beta0)
  names(beta_hat) = 
    c(glue("X_orig{1:ncol(X_orig)}"), glue("X_aug{1:ncol(X_aug)}"))
  
  return(list(intercept_hat = intercept_hat,
              beta_hat = beta_hat, 
              theta_hat_no_intercept = theta_hat_no_intercept,
              converge = NA, 
              iter = NA,
              singular_hessian = NA,
              final_diff = NA,
              phi = mean(mod$phi)))
}


