## ----setup, include = TRUE, echo = FALSE-----------------------
knitr::opts_chunk$set(echo = TRUE, eval = FALSE)
klippy::klippy(position = c("top", "right"), tooltip_message = "Copy")


## ---- message = F, warning = F---------------------------------
library(MASS)
library(tidyverse)
library(adaptBayes)
library(rstan)
library(rstanarm)
library(mnormt)
library(glue)
library(mice)
library(stringr)
library(ggthemes)
library(logistf)


## --------------------------------------------------------------
my_computer = FALSE  
which_batch = 1
convergence_tol = 1e-4
skip_bayes = FALSE 


## --------------------------------------------------------------
if(my_computer) {
  #Choose from between 1-999
  array_id = 25;# 7;#
  # 'sims_per_scenario' is the desired number of independent simulated datasets
  # to be generated. 'total_num_arrays_by_batch' is the total number of SLURM
  # arrays, i.e. array_ids, that will be run for each batch. The code will take
  # care of allocating the appropriate number of arrays.
  sims_per_scenario_by_batch = c(200, 400, 400);
  total_num_arrays_by_batch = c(999, 999, 999); 
} else {
  array_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'));
  sims_per_scenario_by_batch = c(200, 400, 400);
  total_num_arrays_by_batch = c(999, 999, 999); 
}
stopifnot(which_batch <= length(total_num_arrays_by_batch))
stopifnot(length(sims_per_scenario_by_batch) == length(total_num_arrays_by_batch))
sims_per_scenario = sims_per_scenario_by_batch[which_batch]

#Helpful to print warnings when they occur for debugging
options(warn = 1);
#Recommended options from rstan:
options(mc.cores = parallel::detectCores());



## --------------------------------------------------------------
n_hist_seq = c(100, 400, 1600);
n_curr_seq = c(100, 400);
different_external_model_seq = c(F, T);
different_covariate_dist_seq = c(F, T);

source("sim_functions/OR_from_AUC.R")
source("sim_functions/generate_params.R")
source("sim_functions/scenario_setup.R")
source("sim_functions/score_method.R")
source("aux_functions/log1plex.R")
source("aux_functions/get_num_eff.R")


## --------------------------------------------------------------

curr_row <- 
  all_scenarios %>%
  filter(array_id >= start_array_id, 
         array_id <= end_array_id) 
curr_scenario_id <- 
  curr_row %>% pull(scenario)


## --------------------------------------------------------------
if(which_batch == 1) {
  array_id_with_offset = array_id;
} else {
  array_id_with_offset = array_id + sum(total_num_arrays_by_batch[1:(which_batch-1)])
}


## --------------------------------------------------------------
source("methods/glm_vanilla.R");
source("methods/glm_bayes.R");
source("methods/cml_logistic_newtonraphson.R");
source("methods/cml_logistic_saddlepoint.R");
source("methods/gim_logistic_saddlepoint.R");
source("methods/ratios.R")


## --------------------------------------------------------------
which_seeds <- 
  rep(seq(from = all_scenarios %>% slice(curr_scenario_id) %>% pull(start_array_id), 
          to = all_scenarios %>% slice(curr_scenario_id) %>% pull(end_array_id), 
          by = 1), 
      length = sims_per_scenario) == array_id
n_sim = sum(which_seeds)
set.seed(which_batch);
data_seeds <- 
  sample(.Machine$integer.max, sims_per_scenario)[which_seeds]


## --------------------------------------------------------------
true_betas <- c(scenario_list[[pull(curr_row,"scenario_name")]][["orig"]],
                scenario_list[[pull(curr_row,"scenario_name")]][["aug"]])
num_all_coef = length(true_betas)

n_curr <- pull(curr_row, n_curr)
n_hist <- pull(curr_row, n_hist)
different_external_model <- pull(curr_row, different_external_model)
different_covariate_dist <- pull(curr_row, different_covariate_dist)
true_bivariate_cor <- scenario_list[[pull(curr_row,"scenario_name")]]$true_bivariate_cor
true_sigma_x <- true_bivariate_cor + diag(1 - true_bivariate_cor, num_all_coef)
true_chol_sigma_x <- chol(true_sigma_x)

num_orig = length(scenario_list[[pull(curr_row,"scenario_name")]][["orig"]])
num_aug = num_all_coef - num_orig


## --------------------------------------------------------------
set.seed(1)
beta_scale_varying_phi = 
  solve_for_hiershrink_scale(target_mean = get_num_eff(num_all_coef),
                             npar = num_all_coef,
                             local_dof = 1, 
                             global_dof = 1,
                             slab_dof = 4,
                             slab_scale = 2.5,
                             n = n_curr,
                             n_sim = 1e6)$scale %>%
  signif(digits = 4)


## --------------------------------------------------------------
i = 1
all_scores = all_coef_ests = NULL
array_id_stats = 
  curr_row %>% 
  mutate(array_id = array_id,
         n_sim_this_array_id = n_sim, 
         total_runtime_secs = NA)
begin = Sys.time()


## --------------------------------------------------------------
for(i in 1:n_sim) {
  
  
  if(different_external_model) {
    set.seed(data_seeds[i])
    true_betas_hist <- 
      true_betas + rnorm(length(true_betas), sd = mean(abs(true_betas)) / 3)
    true_betas_hist <- 
      round(true_betas_hist * (sum(abs(true_betas)) / sum(abs(true_betas_hist))), 3)
    
  } else {
    true_betas_hist <- true_betas
  }
  
  if(different_covariate_dist) {
    set.seed(data_seeds[i])
    foo <- rWishart(1, 2*num_all_coef^2, true_sigma_x / (2*num_all_coef^2))[,,1]
    foo <- diag(1/sqrt(diag(foo))) %*% foo %*% diag(1/sqrt(diag(foo)))
    true_chol_sigma_x_hist <- chol(foo)
    rm(foo)
  } else {
    true_chol_sigma_x_hist <- true_chol_sigma_x
  }
  
  set.seed(data_seeds[i])
  source("sim_functions/draw_single_dataset.R")
  cat(glue("scenario = {curr_scenario_id}, array_id = {array_id_with_offset}; i/n_sim = {i}/{n_sim}; seed = {data_seeds[i]}\n\n"))
  
  # Thetas (historical study) ----
  
  all_coef_ests <- 
    bind_rows(all_coef_ests, 
              bind_cols(data_seed = data_seeds[i], 
                        sim_num = i, 
                        beta_label = 1:num_all_coef,
                        true_betas = true_betas, 
                        est_betas = c(theta_tilde, rep(NA, num_aug)),
                        method_name = "historical",
                        orig = unlist(map2(c(T,F), c(num_orig, num_aug), rep))))
  
  
  # True Betas ----
  fit_truth <- list(intercept_hat = true_beta0, 
                    beta_hat = true_betas, 
                    iter = NA, 
                    converge = NA, 
                    singular_hessian = NA, 
                    final_diff = NA,
                    objective_function = c("total" = NA))
  all_scores <- 
    bind_rows(all_scores, 
              bind_cols(data_seed =  data_seeds[i], 
                        sim_num = i, 
                        method_name = "truth", 
                        score_method(fit_truth, true_betas, num_orig, num_aug, Xv, Yv)))
  
  rm(fit_truth)
  
  
  # GLM Vanilla ----
  fit_glm_vanilla <- fxn_glm_vanilla(Y = Yc, X_orig = Xc_orig, X_aug = Xc_aug)
  score_glm_vanilla <- score_method(fit_glm_vanilla, true_betas, num_orig, num_aug, Xv, Yv)
  all_scores <- 
    bind_rows(all_scores, 
              bind_cols(data_seed =  data_seeds[i], 
                        sim_num = i, 
                        method_name = "glm_vanilla", 
                        score_glm_vanilla))
  all_coef_ests <- 
    bind_rows(all_coef_ests, 
              bind_cols(data_seed =  data_seeds[i], 
                        sim_num = i, 
                        method_name = "glm_vanilla", 
                        score_glm_vanilla,
                        beta_label = 1:num_all_coef,
                        true_betas = true_betas, 
                        est_betas = fit_glm_vanilla$beta_hat, 
                        orig = unlist(map2(c(T,F), c(num_orig, num_aug), rep))))
  rm(fit_glm_vanilla, score_glm_vanilla)
  
  # GLM Bayes ----
  if(!skip_bayes) {
    fit_glm_bayes <- fxn_glm_bayes(Y = Yc, X_orig = Xc_orig, X_aug = Xc_aug, seed = data_seeds[i])
    score_glm_bayes <- score_method(fit_glm_bayes, true_betas, num_orig, num_aug, Xv, Yv)
    all_scores <- 
      bind_rows(all_scores, 
                bind_cols(data_seed =  data_seeds[i], 
                          sim_num = i, 
                          method_name = "glm_bayes", 
                          score_glm_bayes))
    
    all_coef_ests <- 
      bind_rows(all_coef_ests, 
                bind_cols(data_seed =  data_seeds[i], 
                          sim_num = i, 
                          method_name = "glm_bayes", 
                          score_glm_bayes,
                          beta_label = 1:num_all_coef,
                          true_betas = true_betas, 
                          est_betas = fit_glm_bayes$beta_hat, 
                          orig = unlist(map2(c(T,F), c(num_orig, num_aug), rep))))
    rm(fit_glm_bayes, score_glm_bayes)
  }
  
  # Ratios ----
  fit_ratios <- ratios(Y = Yc, X_orig = Xc_orig, X_aug = Xc_aug, 
                       theta_tilde = theta_tilde)
  score_ratios <- score_method(fit_ratios, true_betas, num_orig, num_aug, Xv, Yv)
  
  all_scores <- 
    bind_rows(all_scores, 
              bind_cols(data_seed =  data_seeds[i], 
                        sim_num = i, 
                        method_name = "ratios", 
                        score_ratios))
  all_coef_ests <- 
    bind_rows(all_coef_ests, 
              bind_cols(data_seed =  data_seeds[i], 
                        sim_num = i, 
                        method_name = "ratios", 
                        score_ratios, 
                        beta_label = 1:num_all_coef,
                        true_betas = true_betas, 
                        est_betas = fit_ratios$beta_hat, 
                        orig = unlist(map2(c(T,F), c(num_orig, num_aug), rep))))
  rm(fit_ratios, score_ratios)
  
  # CML (Saddlepoint) ----
  
  fit_cml_logistic_saddlepoint <- fxn_cml_logistic_saddlepoint(Y = Yc, X_orig = Xc_orig, X_aug = Xc_aug, 
                                                               theta_tilde_with_intercept = theta_tilde_with_intercept,
                                                               beta_init_with_intercept = beta_init_with_intercept,
                                                               tol = convergence_tol, max_rep = 1000, max_lambda = Inf)
  score_cml_logistic_saddlepoint <- score_method(fit_cml_logistic_saddlepoint, true_betas, num_orig, num_aug, Xv, Yv)
  all_scores <- 
    bind_rows(all_scores, 
              bind_cols(data_seed =  data_seeds[i], 
                        sim_num = i,
                        method_name = "cml_saddlepoint", 
                        score_cml_logistic_saddlepoint))
  all_coef_ests <- 
    bind_rows(all_coef_ests, 
              bind_cols(data_seed =  data_seeds[i], 
                        sim_num = i, 
                        method_name = "cml_saddlepoint", 
                        score_cml_logistic_saddlepoint,
                        beta_label = 1:num_all_coef,
                        true_betas = true_betas, 
                        est_betas = fit_cml_logistic_saddlepoint$beta_hat, 
                        orig = unlist(map2(c(T,F), c(num_orig, num_aug), rep))))
  rm(fit_cml_logistic_saddlepoint, score_cml_logistic_saddlepoint)
  
  
  # CML (Newton Raphson) ----
  fit_cml_logistic_newtonraphson <- fxn_cml_logistic_newtonraphson(Y = Yc, X_orig = Xc_orig, X_aug = Xc_aug, 
                                                                   theta_tilde_with_intercept = theta_tilde_with_intercept,
                                                                   beta_init_with_intercept = beta_init_with_intercept, 
                                                                   tol = convergence_tol, max_rep = 1000)
  score_cml_logistic_newtonraphson <- score_method(fit_cml_logistic_newtonraphson, true_betas, num_orig, num_aug, Xv, Yv)
  all_scores <- 
    bind_rows(all_scores, 
              bind_cols(data_seed =  data_seeds[i], 
                        sim_num = i,
                        method_name = "cml_newtonraphson", 
                        score_cml_logistic_newtonraphson))
  all_coef_ests <- 
    bind_rows(all_coef_ests, 
              bind_cols(data_seed =  data_seeds[i], 
                        sim_num = i, 
                        method_name = "cml_newtonraphson", 
                        score_cml_logistic_newtonraphson,
                        beta_label = 1:num_all_coef,
                        true_betas = true_betas, 
                        est_betas = fit_cml_logistic_newtonraphson$beta_hat, 
                        orig = unlist(map2(c(T,F), c(num_orig, num_aug), rep))))
  rm(fit_cml_logistic_newtonraphson, score_cml_logistic_newtonraphson)
  
  # CML (Newton Raphson, stepwise) ----
  fit_cml_logistic_newtonraphson_step <- fxn_cml_logistic_newtonraphson(Y = Yc, X_orig = Xc_orig, X_aug = Xc_aug, 
                                                                        theta_tilde_with_intercept = theta_tilde_with_intercept,
                                                                        beta_init_with_intercept = beta_init_with_intercept, 
                                                                        tol = convergence_tol, max_rep = 1000, step = 0.1)
  score_cml_logistic_newtonraphson_step <- score_method(fit_cml_logistic_newtonraphson_step, true_betas, num_orig, num_aug, Xv, Yv)
  all_scores <- 
    bind_rows(all_scores, 
              bind_cols(data_seed =  data_seeds[i], 
                        sim_num = i,
                        method_name = "cml_newtonraphson_step", 
                        score_cml_logistic_newtonraphson_step))
  all_coef_ests <- 
    bind_rows(all_coef_ests, 
              bind_cols(data_seed =  data_seeds[i], 
                        sim_num = i, 
                        method_name = "cml_newtonraphson_step",
                        score_cml_logistic_newtonraphson_step,
                        beta_label = 1:num_all_coef,
                        true_betas = true_betas, 
                        est_betas = fit_cml_logistic_newtonraphson_step$beta_hat, 
                        orig = unlist(map2(c(T,F), c(num_orig, num_aug), rep))))
  rm(fit_cml_logistic_newtonraphson_step, score_cml_logistic_newtonraphson_step)
  
  # GIM (original proposal) ----
  
  fit_gim_logistic_saddlepoint_orig <- fxn_gim_logistic_saddlepoint(Y = Yc, X_orig = Xc_orig, X_aug = Xc_aug, 
                                                                    theta_tilde_with_intercept = theta_tilde_with_intercept, 
                                                                    beta_init_with_intercept = beta_init_with_intercept, 
                                                                    n_hist = n_hist, 
                                                                    Sigma0 = NULL,
                                                                    tol = convergence_tol, max_rep = 1000, max_lambda = Inf);
  score_gim_logistic_saddlepoint_orig <- score_method(fit_gim_logistic_saddlepoint_orig, true_betas, num_orig, num_aug, Xv, Yv)
  
  all_scores <- 
    bind_rows(all_scores, 
              bind_cols(data_seed =  data_seeds[i], 
                        sim_num = i, 
                        method_name = "gim_orig", 
                        score_gim_logistic_saddlepoint_orig))
  all_coef_ests <- 
    bind_rows(all_coef_ests, 
              bind_cols(data_seed =  data_seeds[i], 
                        sim_num = i, 
                        method_name = "gim_orig", 
                        score_gim_logistic_saddlepoint_orig,
                        beta_label = 1:num_all_coef,
                        true_betas = true_betas, 
                        est_betas = fit_gim_logistic_saddlepoint_orig$beta_hat, 
                        orig = unlist(map2(c(T,F), c(num_orig, num_aug), rep))))
  rm(fit_gim_logistic_saddlepoint_orig, score_gim_logistic_saddlepoint_orig)
  
  # SAB, historical value of var_theta_tilde ----
  # beta^o + beta^a ~ N(theta_tilde, eta * var_theta_tilde)
  if(!skip_bayes) {
    eigendecomp_hist_var = eigen(var_theta_tilde);
    aug_projection = create_projection(x_curr_orig = Xc_orig %>% `colnames<-`(as.character(glue("p{1:num_orig}"))),
                                       x_curr_aug = Xc_aug %>% `colnames<-`(as.character(glue("q{1:num_aug}"))),
                                       eigenvec_hist_var = t(eigendecomp_hist_var$vectors),
                                       imputes_list = list(c(1, 500)))[[1]];
    sab <-
      glm_sab(y = Yc, 
              x_standardized = cbind(Xc_orig, Xc_aug), 
              family = "binomial",
              alpha_prior_mean = theta_tilde,
              alpha_prior_cov = var_theta_tilde,
              aug_projection = aug_projection,
              phi_dist = "trunc_norm",
              phi_mean = 1,
              phi_sd = 0.25,
              eta_param = 2.5,
              beta_orig_scale = beta_scale_varying_phi, 
              beta_aug_scale = beta_scale_varying_phi, 
              local_dof = 1, 
              global_dof = 1, 
              slab_dof = 4, 
              slab_scale = 2.5,
              only_prior = F, 
              mc_warmup = 1e3, 
              mc_iter_after_warmup = 1e3, 
              mc_chains = 2, 
              mc_thin = 1,
              mc_stepsize = 0.1,
              mc_adapt_delta = 0.999,
              mc_max_treedepth = 15,
              eigendecomp_hist_var = eigendecomp_hist_var,
              seed = data_seeds[i]);
    
    fit_sab <-
      list(beta_hat = apply(sab$beta, 2, median) %>% `names<-`(c(glue("Xc_orig{1:num_orig}"), glue("Xc_aug{1:num_aug}"))),
           intercept_hat = median(sab$mu),
           converge = NA, 
           iter = NA,
           singular_hessian = NA,
           final_diff = NA,
           objective_function = c("total" = NA),
           phi = median(sab$phi), 
           eta = median(sab$eta))
    score_sab <- score_method(fit_sab, true_betas, num_orig, num_aug, Xv, Yv)
    
    all_scores <- 
      bind_rows(all_scores, 
                bind_cols(data_seed =  data_seeds[i], 
                          sim_num = i, 
                          method_name = "sab", 
                          score_sab))
    all_coef_ests <- 
      bind_rows(all_coef_ests, 
                bind_cols(data_seed =  data_seeds[i], 
                          sim_num = i, 
                          method_name = "sab", 
                          score_sab, 
                          beta_label = 1:num_all_coef,
                          true_betas = true_betas, 
                          est_betas = fit_sab$beta_hat, 
                          orig = unlist(map2(c(T,F), c(num_orig, num_aug), rep))))
    rm(sab, fit_sab, score_sab)
    rm(aug_projection)
    
  }
  # SAB2, lm_proj, omega_sq_in_variance = TRUE ----
  # beta^o + beta^a ~ N(omega * theta_tilde, omega^2 * var_theta_tilde)
  
  if(!skip_bayes) {
    
    # Estimate the X_aug ~ X_orig multivariate regression coefficients using
    # weakly penalized ridge regression: beta ~ N(0, 25 * sigma^2)
    lambda_matrix = diag(c(0, rep(0.04, num_orig)))
    aug_projection = as.matrix(coef(lm(rbind(Xc_aug,matrix(0, nrow = num_orig + 1, ncol = num_aug)) ~ -1 + rbind(cbind(1, Xc_orig), lambda_matrix))))[-1,,drop = FALSE]
    dimnames(aug_projection) = NULL;
    rm(lambda_matrix)
    
    sab2 <-
      glm_sab2(y = Yc, 
               x_standardized = cbind(Xc_orig, Xc_aug), 
               family = "binomial",
               alpha_prior_mean = theta_tilde,
               alpha_prior_cov = var_theta_tilde,
               aug_projection = aug_projection,
               phi_dist = "trunc_norm",
               phi_mean = 1,
               phi_sd = 0.25,
               eta_param = Inf,
               omega_mean = 0, 
               omega_sd = 0.25,
               omega_sq_in_variance = TRUE,
               beta_orig_scale = beta_scale_varying_phi, 
               beta_aug_scale = beta_scale_varying_phi, 
               local_dof = 1, 
               global_dof = 1, 
               slab_dof = 4, 
               slab_scale = 2.5,
               only_prior = F, 
               mc_warmup = 1e3, 
               mc_iter_after_warmup = 1e3, 
               mc_chains = 2, 
               mc_thin = 1,
               mc_stepsize = 0.1,
               mc_adapt_delta = 0.999,
               mc_max_treedepth = 15,
               eigendecomp_hist_var = eigendecomp_hist_var,
               seed = data_seeds[i]);
    
    fit_sab2 <-
      list(beta_hat = apply(sab2$beta, 2, median) %>% `names<-`(c(glue("Xc_orig{1:num_orig}"), glue("Xc_aug{1:num_aug}"))),
           intercept_hat = median(sab2$mu),
           converge = NA, 
           iter = NA,
           singular_hessian = NA,
           final_diff = NA,
           objective_function = c("total" = NA),
           phi = median(sab2$phi), 
           eta = median(sab2$eta),
           omega = median(sab2$omega))
    score_sab2 <- score_method(fit_sab2, true_betas, num_orig, num_aug, Xv, Yv)
    
    all_scores <- 
      bind_rows(all_scores, 
                bind_cols(data_seed =  data_seeds[i], 
                          sim_num = i, 
                          method_name = "sab2", 
                          score_sab2))
    all_coef_ests <- 
      bind_rows(all_coef_ests, 
                bind_cols(data_seed = data_seeds[i], 
                          sim_num = i, 
                          method_name = "sab2",
                          score_sab2,
                          beta_label = 1:num_all_coef,
                          true_betas = true_betas, 
                          est_betas = fit_sab2$beta_hat, 
                          orig = unlist(map2(c(T,F), c(num_orig, num_aug), rep))))
    rm(sab2, fit_sab2, score_sab2)
  }
  
  # SAB2, lm_proj, omega_sq_in_variance = FALSE ----
  # beta^o + beta^a ~ N(omega * theta_tilde, eta * var_theta_tilde)
  if(!skip_bayes) {
    sab2 <-
      glm_sab2(y = Yc, 
               x_standardized = cbind(Xc_orig, Xc_aug), 
               family = "binomial",
               alpha_prior_mean = theta_tilde,
               alpha_prior_cov = var_theta_tilde,
               aug_projection = aug_projection,
               phi_dist = "trunc_norm",
               phi_mean = 1,
               phi_sd = 0.25,
               eta_param = 2.5, # 1st distinguishing feature from above version of SAB2
               omega_mean = 0, 
               omega_sd = 0.25,
               omega_sq_in_variance = FALSE, # 2st distinguishing feature from above version of SAB2
               beta_orig_scale = beta_scale_varying_phi, 
               beta_aug_scale = beta_scale_varying_phi, 
               local_dof = 1, 
               global_dof = 1, 
               slab_dof = 4, 
               slab_scale = 2.5,
               only_prior = F, 
               mc_warmup = 1e3, 
               mc_iter_after_warmup = 1e3, 
               mc_chains = 2, 
               mc_thin = 1,
               mc_stepsize = 0.1,
               mc_adapt_delta = 0.999,
               mc_max_treedepth = 15,
               eigendecomp_hist_var = eigendecomp_hist_var,
               seed = data_seeds[i]);
    
    fit_sab2 <-
      list(beta_hat = apply(sab2$beta,2,median) %>% `names<-`(c(glue("Xc_orig{1:num_orig}"), glue("Xc_aug{1:num_aug}"))),
           intercept_hat = median(sab2$mu),
           converge = NA, 
           iter = NA,
           singular_hessian = NA,
           final_diff = NA,
           objective_function = c("total" = NA),
           phi = median(sab2$phi), 
           eta = median(sab2$eta),
           omega = median(sab2$omega))
    score_sab2 <- score_method(fit_sab2, true_betas, num_orig, num_aug, Xv, Yv)
    
    all_scores <- 
      bind_rows(all_scores, 
                bind_cols(data_seed = data_seeds[i], 
                          sim_num = i,
                          method_name = "sab2_alt", 
                          score_sab2))
    all_coef_ests <- 
      bind_rows(all_coef_ests, 
                bind_cols(data_seed = data_seeds[i], 
                          sim_num = i, 
                          method_name = "sab2_alt",
                          score_sab2,
                          beta_label = 1:num_all_coef,
                          true_betas = true_betas, 
                          est_betas = fit_sab2$beta_hat, 
                          orig = unlist(map2(c(T,F), c(num_orig, num_aug), rep))))
    rm(sab2, fit_sab2, score_sab2)
  }
  rm(var_theta_tilde);
  rm(eigendecomp_hist_var)
  
  # SAB2, lm_proj, omega_sq_in_variance = TRUE, internal estimate of var_theta ----
  # beta^o + beta^a ~ N(omega * theta_tilde, omega^2 * var_theta_tilde)
  if(!skip_bayes) {
    
    var_theta_tilde_alt <- (n_curr / n_hist) * vcov(logistf(Yh ~ Xh_orig, family="binomial", plconf = 1))[-1, -1, drop = FALSE]
    
    sab2 <-
      glm_sab2(y = Yc, 
               x_standardized = cbind(Xc_orig, Xc_aug), 
               family = "binomial",
               alpha_prior_mean = theta_tilde,
               alpha_prior_cov = var_theta_tilde_alt,
               aug_projection = aug_projection,
               phi_dist = "trunc_norm",
               phi_mean = 1,
               phi_sd = 0.25,
               eta_param = Inf,
               omega_mean = 0, 
               omega_sd = 0.25,
               omega_sq_in_variance = TRUE,
               beta_orig_scale = beta_scale_varying_phi, 
               beta_aug_scale = beta_scale_varying_phi, 
               local_dof = 1, 
               global_dof = 1, 
               slab_dof = 4, 
               slab_scale = 2.5,
               only_prior = F, 
               mc_warmup = 1e3, 
               mc_iter_after_warmup = 1e3, 
               mc_chains = 2, 
               mc_thin = 1,
               mc_stepsize = 0.1,
               mc_adapt_delta = 0.999,
               mc_max_treedepth = 15,
               eigendecomp_hist_var = NULL,
               seed = data_seeds[i]);
    
    fit_sab2 <-
      list(beta_hat = apply(sab2$beta,2,median) %>% `names<-`(c(glue("Xc_orig{1:num_orig}"), glue("Xc_aug{1:num_aug}"))),
           intercept_hat = median(sab2$mu),
           converge = NA, 
           iter = NA,
           singular_hessian = NA,
           final_diff = NA,
           objective_function = c("total" = NA),
           phi = median(sab2$phi), 
           eta = median(sab2$eta),
           omega = median(sab2$omega))
    score_sab2 <- score_method(fit_sab2, true_betas, num_orig, num_aug, Xv, Yv)
    
    all_scores <- 
      bind_rows(all_scores, 
                bind_cols(data_seed =  data_seeds[i],
                          sim_num = i,
                          method_name = "sab2_flexible", 
                          score_sab2))
    all_coef_ests <- 
      bind_rows(all_coef_ests, 
                bind_cols(data_seed =  data_seeds[i], 
                          sim_num = i, 
                          method_name = "sab2_flexible",
                          score_sab2,
                          beta_label = 1:num_all_coef,
                          true_betas = true_betas, 
                          est_betas = fit_sab2$beta_hat, 
                          orig = unlist(map2(c(T,F), c(num_orig, num_aug), rep))))
    rm(sab2, fit_sab2, score_sab2)
    rm(aug_projection)
    rm(var_theta_tilde_alt)
  }
  
  array_id_stats <-
    array_id_stats %>%
    mutate(total_runtime_secs = as.numeric(difftime(Sys.time(), begin, units =  "secs")))
  
  rm(beta_init_with_intercept, mod_hist, 
     true_betas_hist,
     true_chol_sigma_x_hist,
     Pc, Ph, 
     theta_tilde_with_intercept, 
     theta_tilde, 
     var_theta_tilde_with_intercept, 
     Xc, Xc_aug, Xc_orig, 
     Xh, Xh_aug, Xh_orig,     
     Yc, Yh, 
     Xv, Pv, Yv)
  
  all_scores <- 
    all_scores %>%
    mutate(array_id = array_id_with_offset,
           which_batch = which_batch, 
           n_hist = n_hist, 
           n_curr = n_curr, 
           different_external_model = different_external_model,
           different_covariate_dist = different_covariate_dist,
           true_bivariate_cor = true_bivariate_cor, 
           scenario_name = pull(curr_row,"scenario_name")) %>%
    fill(array_id, which_batch, n_hist, n_curr, true_bivariate_cor, 
         different_external_model, different_covariate_dist) %>%
    select(array_id, sim_num, which_batch, data_seed,
           n_hist, n_curr, true_bivariate_cor, 
           different_external_model, different_covariate_dist, 
           scenario_name, 
           everything())
  
  all_coef_ests <- 
    all_coef_ests %>%
    mutate(array_id = array_id_with_offset,
           which_batch = which_batch, 
           n_hist = n_hist, 
           n_curr = n_curr, 
           different_external_model = different_external_model,
           different_covariate_dist = different_covariate_dist,
           true_bivariate_cor = true_bivariate_cor, 
           scenario_name = pull(curr_row,"scenario_name")) %>%
    fill(array_id, which_batch, n_hist, n_curr, true_bivariate_cor, 
         different_external_model, different_covariate_dist) %>%
    select(array_id, sim_num, which_batch, data_seed,
           n_hist, n_curr, true_bivariate_cor, 
           different_external_model, different_covariate_dist, 
           scenario_name, 
           method_name, true_betas, orig,
           everything())
  
  
  if(!my_computer) {
    # Save scores and estimates and running times as csv
    write_csv(all_scores,
              file = paste0("out/job",array_id_with_offset, "_scores.csv"),
              append = FALSE);
    
    write_csv(all_coef_ests,
              file = paste0("out/job",array_id_with_offset, "_coefs.csv"),
              append = FALSE);
    
    write_csv(array_id_stats,
              file = paste0("out/job",array_id_with_offset, "_array_stats.csv"),
              append = FALSE);
    
  }
}

if(my_computer) {
  all_coef_ests %>% 
    filter(beta_label == 1) %>%
    select(sim_num, method_name, est_betas) %>% 
    pivot_wider(id_cols = sim_num, names_from = method_name, values_from = est_betas)
  all_scores %>% 
    group_by(method_name) %>% 
    summarize(mean(squared_error), mean(squared_error_orig), mean(squared_error_aug))
}


