## ----message = F, warning = F---------------------------------
library(MASS)
library(tidyverse)
library(adaptBayes)
library(rstan)
library(rstanarm)
library(glue)
library(mice)
library(stringr)
library(ggthemes)
library(logistf)
library(GENMETA)
library(gim)


## -------------------------------------------------------------
my_computer = FALSE
which_batch = 3
convergence_tol = 1e-4
skip_bayes = FALSE


## -------------------------------------------------------------
if(my_computer) {
  #Choose from between 1-999
  array_id = 1;# 7;#data_seed=1140350788
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



## -------------------------------------------------------------
n_hist_seq = c(400, 1600);
n_curr_seq = c(200, 400);
different_external_model_seq = c(F, T);
different_covariate_dist_seq = c(F, T);

source("sim_functions/OR_from_AUC.R")
source("sim_functions/generate_params.R")
source("sim_functions/scenario_setup.R")
source("sim_functions/score_method.R")
source("aux_functions/log1plex.R")
source("aux_functions/get_num_eff.R")


## -------------------------------------------------------------

curr_row <- 
  all_scenarios %>%
  filter(array_id >= start_array_id, 
         array_id <= end_array_id) 
curr_scenario_id <- 
  curr_row %>% pull(scenario)


## -------------------------------------------------------------
if(which_batch == 1) {
  array_id_with_offset = array_id;
} else {
  array_id_with_offset = array_id + sum(total_num_arrays_by_batch[1:(which_batch-1)])
}


## -------------------------------------------------------------
source("methods/glm_vanilla.R");
source("methods/glm_bayes.R");
source("methods/cml_logistic_newtonraphson.R");
source("methods/cml_logistic_saddlepoint.R");
source("methods/gim_logistic_saddlepoint.R");
source("methods/gim_logistic_author.R");
source("methods/genmeta_logistic_author.R");
source("methods/ratios.R")


## -------------------------------------------------------------
which_seeds <- 
  rep(seq(from = all_scenarios %>% slice(curr_scenario_id) %>% pull(start_array_id), 
          to = all_scenarios %>% slice(curr_scenario_id) %>% pull(end_array_id), 
          by = 1), 
      length = sims_per_scenario) == array_id
n_sim = sum(which_seeds)
set.seed(which_batch);
data_seeds <- 
  sample(.Machine$integer.max, sims_per_scenario)[which_seeds]


## -------------------------------------------------------------
num_orig = length(scenario_list[[pull(curr_row,"scenario_name")]][["orig"]])
orig_var_names = as.character(glue("p{1:num_orig}"))
num_aug = length(scenario_list[[pull(curr_row,"scenario_name")]][["aug"]])
aug_var_names = as.character(glue("q{1:num_aug}"))
num_all = num_orig + num_aug;
all_var_names = c(orig_var_names, aug_var_names);

n_curr <- pull(curr_row, n_curr)
n_hist <- pull(curr_row, n_hist)
different_external_model <- pull(curr_row, different_external_model)
different_covariate_dist <- pull(curr_row, different_covariate_dist)
true_bivariate_cor <- scenario_list[[pull(curr_row,"scenario_name")]][["true_bivariate_cor"]]
true_sigma_x <- true_bivariate_cor + diag(1 - true_bivariate_cor, num_all)
true_chol_sigma_x <- chol(true_sigma_x)
true_betas <- 
  c(scenario_list[[pull(curr_row,"scenario_name")]][["orig"]],
    scenario_list[[pull(curr_row,"scenario_name")]][["aug"]]) %>%
  `names<-`(all_var_names)
which_binary <- scenario_list[[pull(curr_row,"scenario_name")]][["which_binary"]]


## -------------------------------------------------------------
slab_dof = 4;
slab_scale = 2.5;
set.seed(1)
beta_scale_varying_phi = 
  solve_for_hiershrink_scale(target_mean = get_num_eff(num_all),
                             npar = num_all,
                             local_dof = 1, 
                             global_dof = 1,
                             slab_dof = slab_dof,
                             slab_scale = slab_scale,
                             n = n_curr,
                             n_sim = 1e6)$scale %>%
  signif(digits = 4)


## -------------------------------------------------------------
i = 1
all_scores = all_beta_ests = all_lambda_ests = all_bayesian_diag = NULL
array_id_stats = 
  curr_row %>% 
  mutate(array_id = array_id_with_offset,
         n_sim_this_array_id = n_sim, 
         total_runtime_secs = NA)


## -------------------------------------------------------------
begin = Sys.time()
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
    foo <- rWishart(1, 2*num_all^2, true_sigma_x / (2*num_all^2))[,,1]
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
  all_beta_ests <- 
    bind_rows(all_beta_ests, 
              bind_cols(data_seed = data_seeds[i], 
                        sim_num = i, 
                        beta_label = all_var_names,
                        true_betas = true_betas, 
                        est_betas = c(theta_tilde, rep(NA, num_aug)) %>% `names<-`(all_var_names),
                        method_name = "historical",
                        orig = unlist(map2(c(T,F), c(num_orig, num_aug), rep))))
  
  
  # True Betas ----
  begin2 <- Sys.time();
  curr_fit <- list(intercept_hat = true_beta0, 
                   beta_hat = true_betas, 
                   converge = NA, 
                   message = NA, 
                   iter = NA, 
                   singular_hessian = NA, 
                   final_diff = NA,
                   objective_function = c("total" = NA))
  curr_run_time <- as.numeric(difftime(Sys.time(), begin2, units =  "secs"));
  
  all_scores <- 
    bind_rows(all_scores, 
              bind_cols(data_seed =  data_seeds[i], 
                        sim_num = i, 
                        method_name = "truth", 
                        run_time = curr_run_time,
                        score_method(curr_fit, true_betas, num_orig, num_aug, Xv, Yv)))
  
  rm(curr_fit, curr_run_time)
  
  
  # GLM Vanilla ----
  begin2 <- Sys.time();
  curr_fit <- fxn_glm_vanilla(Y = Yc, X = Xc)
  curr_run_time <- as.numeric(difftime(Sys.time(), begin2, units =  "secs"));
  
  curr_score <- score_method(curr_fit, true_betas, num_orig, num_aug, Xv, Yv)
  all_scores <- 
    bind_rows(all_scores, 
              bind_cols(data_seed =  data_seeds[i], 
                        sim_num = i, 
                        method_name = "glm_vanilla", 
                        run_time = curr_run_time,
                        curr_score))
  all_beta_ests <- 
    bind_rows(all_beta_ests, 
              bind_cols(data_seed =  data_seeds[i], 
                        sim_num = i, 
                        method_name = "glm_vanilla", 
                        curr_score,
                        beta_label = all_var_names,
                        true_betas = true_betas, 
                        est_betas = curr_fit$beta_hat[all_var_names], 
                        orig = unlist(map2(c(T,F), c(num_orig, num_aug), rep))))
  # Keep this for the Generalized Meta Analysis
  standard_mle <- 
    c(curr_fit$intercept_hat, curr_fit$beta_hat[all_var_names])
  rm(curr_fit, curr_run_time, curr_score)
  
  # GLM Bayes ----
  if(!skip_bayes) {
    begin2 <- Sys.time();
    curr_fit <- fxn_glm_bayes(Y = Yc, X = Xc, seed = data_seeds[i])
    curr_run_time <- as.numeric(difftime(Sys.time(), begin2, units =  "secs"));
    curr_score <- score_method(curr_fit, true_betas, num_orig, num_aug, Xv, Yv)
    all_scores <- 
      bind_rows(all_scores, 
                bind_cols(data_seed = data_seeds[i], 
                          sim_num = i, 
                          method_name = "glm_bayes", 
                          run_time = curr_run_time,
                          curr_score))
    
    all_beta_ests <- 
      bind_rows(all_beta_ests, 
                bind_cols(data_seed =  data_seeds[i], 
                          sim_num = i, 
                          method_name = "glm_bayes", 
                          curr_score,
                          beta_label = all_var_names,
                          true_betas = true_betas, 
                          est_betas = curr_fit$beta_hat[all_var_names], 
                          orig = unlist(map2(c(T,F), c(num_orig, num_aug), rep))))
    
    all_bayesian_diag <- 
      bind_rows(all_bayesian_diag, 
                bind_cols(data_seed =  data_seeds[i], 
                          sim_num = i, 
                          method_name = "glm_bayes", 
                          num_divergences = curr_fit$num_divergences,
                          num_max_treedepth = NA_real_, # didn't bother to calculate this for this method
                          min_ebfmi = NA_real_, # didn't bother to calculate this for this method
                          max_rhat = curr_fit$max_rhat))
    # Keep this for the Generalized Meta Analysis
    #standard_bayes <- 
    #  c(curr_fit$intercept_hat, curr_fit$beta_hat[all_var_names])
    rm(curr_fit, curr_run_time, curr_score)
  }
  
  # Ratios ----
  begin2 <- Sys.time();
  curr_fit <- ratios(Y = Yc, X = Xc, theta_tilde = theta_tilde)
  curr_run_time <- as.numeric(difftime(Sys.time(), begin2, units =  "secs"));
  curr_score <- score_method(curr_fit, true_betas, num_orig, num_aug, Xv, Yv)
  
  all_scores <- 
    bind_rows(all_scores, 
              bind_cols(data_seed =  data_seeds[i], 
                        sim_num = i, 
                        method_name = "ratios", 
                        run_time = curr_run_time,
                        curr_score))
  all_beta_ests <- 
    bind_rows(all_beta_ests, 
              bind_cols(data_seed =  data_seeds[i], 
                        sim_num = i, 
                        method_name = "ratios", 
                        curr_score, 
                        beta_label = all_var_names,
                        true_betas = true_betas, 
                        est_betas = curr_fit$beta_hat[all_var_names], 
                        orig = unlist(map2(c(T,F), c(num_orig, num_aug), rep))))
  rm(curr_fit, curr_run_time, curr_score)
  
  # CML (Saddlepoint) ----
  
  begin2 <- Sys.time();
  curr_fit <- fxn_cml_logistic_saddlepoint(Y = Yc, X = Xc,
                                           theta_tilde_with_intercept = theta_tilde_with_intercept,
                                           beta_init_with_intercept = beta_init_with_intercept,
                                           tol = convergence_tol)
  curr_run_time <- as.numeric(difftime(Sys.time(), begin2, units =  "secs"));
  curr_score <- score_method(curr_fit, true_betas, num_orig, num_aug, Xv, Yv)
  all_scores <- 
    bind_rows(all_scores, 
              bind_cols(data_seed =  data_seeds[i], 
                        sim_num = i,
                        method_name = "cml_saddlepoint", 
                        run_time = curr_run_time,
                        curr_score))
  all_beta_ests <- 
    bind_rows(all_beta_ests, 
              bind_cols(data_seed =  data_seeds[i], 
                        sim_num = i, 
                        method_name = "cml_saddlepoint", 
                        curr_score,
                        beta_label = all_var_names,
                        true_betas = true_betas, 
                        est_betas = curr_fit$beta_hat[all_var_names], 
                        orig = unlist(map2(c(T,F), c(num_orig, num_aug), rep))))
  all_lambda_ests <- 
    bind_rows(all_lambda_ests, 
              bind_cols(data_seed =  data_seeds[i], 
                        sim_num = i, 
                        method_name = "cml_saddlepoint", 
                        lambda_label = c("(Intercept)", orig_var_names),
                        est_lambda = curr_fit$lambda_hat))
  rm(curr_fit, curr_run_time, curr_score)
  
  
  # CML (Newton Raphson) ----
  begin2 <- Sys.time();
  curr_fit <- fxn_cml_logistic_newtonraphson(Y = Yc, X = Xc,
                                             theta_tilde_with_intercept = theta_tilde_with_intercept,
                                             beta_init_with_intercept = beta_init_with_intercept, 
                                             tol = convergence_tol)
  curr_run_time <- as.numeric(difftime(Sys.time(), begin2, units =  "secs"));
  curr_score <- score_method(curr_fit, true_betas, num_orig, num_aug, Xv, Yv)
  all_scores <- 
    bind_rows(all_scores, 
              bind_cols(data_seed =  data_seeds[i], 
                        sim_num = i,
                        method_name = "cml_newtonraphson", 
                        run_time = curr_run_time,
                        curr_score))
  all_beta_ests <- 
    bind_rows(all_beta_ests, 
              bind_cols(data_seed =  data_seeds[i], 
                        sim_num = i, 
                        method_name = "cml_newtonraphson", 
                        curr_score,
                        beta_label = all_var_names,
                        true_betas = true_betas, 
                        est_betas = curr_fit$beta_hat[all_var_names], 
                        orig = unlist(map2(c(T,F), c(num_orig, num_aug), rep))))
  all_lambda_ests <- 
    bind_rows(all_lambda_ests, 
              bind_cols(data_seed =  data_seeds[i], 
                        sim_num = i, 
                        method_name = "cml_newtonraphson", 
                        lambda_label = c("(Intercept)", orig_var_names),
                        est_lambda = curr_fit$lambda_hat))
  rm(curr_fit, curr_run_time, curr_score)
  
  # CML (Newton Raphson, stepwise) ----
  begin2 <- Sys.time();
  curr_fit <- fxn_cml_logistic_newtonraphson(Y = Yc, X = Xc,
                                             theta_tilde_with_intercept = theta_tilde_with_intercept,
                                             beta_init_with_intercept = beta_init_with_intercept, 
                                             tol = convergence_tol, step = 0.1)
  curr_run_time <- as.numeric(difftime(Sys.time(), begin2, units =  "secs"));
  curr_score <- score_method(curr_fit, true_betas, num_orig, num_aug, Xv, Yv)
  all_scores <- 
    bind_rows(all_scores, 
              bind_cols(data_seed =  data_seeds[i], 
                        sim_num = i,
                        method_name = "cml_newtonraphson_step", 
                        run_time = curr_run_time,
                        curr_score))
  all_beta_ests <- 
    bind_rows(all_beta_ests, 
              bind_cols(data_seed =  data_seeds[i], 
                        sim_num = i, 
                        method_name = "cml_newtonraphson_step",
                        curr_score,
                        beta_label = all_var_names,
                        true_betas = true_betas, 
                        est_betas = curr_fit$beta_hat[all_var_names], 
                        orig = unlist(map2(c(T,F), c(num_orig, num_aug), rep))))
  all_lambda_ests <- 
    bind_rows(all_lambda_ests, 
              bind_cols(data_seed =  data_seeds[i], 
                        sim_num = i, 
                        method_name = "cml_newtonraphson_step", 
                        lambda_label = c("(Intercept)", orig_var_names),
                        est_lambda = curr_fit$lambda_hat))
  
  rm(curr_fit, curr_run_time, curr_score)
  
  # GIM (original proposal) ----
  begin2 <- Sys.time();
  curr_fit <- fxn_gim_logistic_saddlepoint(Y = Yc, X = Xc, 
                                           theta_tilde_with_intercept = theta_tilde_with_intercept, 
                                           beta_init_with_intercept = beta_init_with_intercept, 
                                           n_hist = n_hist, 
                                           tol = convergence_tol);
  curr_run_time <- as.numeric(difftime(Sys.time(), begin2, units =  "secs"));
  curr_score <- score_method(curr_fit, true_betas, num_orig, num_aug, Xv, Yv)
  
  all_scores <- 
    bind_rows(all_scores, 
              bind_cols(data_seed = data_seeds[i], 
                        sim_num = i, 
                        method_name = "gim_sandwich", 
                        run_time = curr_run_time,
                        curr_score))
  all_beta_ests <- 
    bind_rows(all_beta_ests, 
              bind_cols(data_seed = data_seeds[i], 
                        sim_num = i, 
                        method_name = "gim_sandwich", 
                        curr_score,
                        beta_label = all_var_names,
                        true_betas = true_betas, 
                        est_betas = curr_fit$beta_hat[all_var_names], 
                        orig = unlist(map2(c(T,F), c(num_orig, num_aug), rep))))
  all_lambda_ests <- 
    bind_rows(all_lambda_ests, 
              bind_cols(data_seed =  data_seeds[i], 
                        sim_num = i, 
                        method_name = "gim_sandwich", 
                        lambda_label = c("(Intercept)", orig_var_names),
                        est_lambda = curr_fit$lambda_hat))
  
  rm(curr_fit, curr_run_time, curr_score)
  
  if(0) {
    # GIM (likelihood-based) ----
    begin2 <- Sys.time();
    curr_fit <- fxn_gim_logistic_saddlepoint(Y = Yc, X = Xc, 
                                             theta_tilde_with_intercept = theta_tilde_with_intercept, 
                                             beta_init_with_intercept = beta_init_with_intercept, 
                                             n_hist = n_hist, 
                                             Sigma_h = "likelihood",
                                             tol = convergence_tol);
    curr_run_time <- as.numeric(difftime(Sys.time(), begin2, units =  "secs"));
    curr_score <- score_method(curr_fit, true_betas, num_orig, num_aug, Xv, Yv)
    
    all_scores <- 
      bind_rows(all_scores, 
                bind_cols(data_seed = data_seeds[i], 
                          sim_num = i, 
                          method_name = "gim_lik", 
                          run_time = curr_run_time,
                          curr_score))
    all_beta_ests <- 
      bind_rows(all_beta_ests, 
                bind_cols(data_seed = data_seeds[i], 
                          sim_num = i, 
                          method_name = "gim_lik", 
                          curr_score,
                          beta_label = all_var_names,
                          true_betas = true_betas, 
                          est_betas = curr_fit$beta_hat[all_var_names], 
                          orig = unlist(map2(c(T,F), c(num_orig, num_aug), rep))))
    all_lambda_ests <- 
      bind_rows(all_lambda_ests, 
                bind_cols(data_seed =  data_seeds[i], 
                          sim_num = i, 
                          method_name = "gim_lik", 
                          lambda_label = c("(Intercept)", orig_var_names),
                          est_lambda = curr_fit$lambda_hat))
    
    rm(curr_fit, curr_run_time, curr_score)
  }
  
  # GIM (using R package gim) ----
  begin2 <- Sys.time();
  curr_fit <- fxn_gim_logistic_author(Y = Yc, X = Xc, 
                                      theta_tilde_with_intercept = theta_tilde_with_intercept, 
                                      n_hist = n_hist);
  curr_run_time <- as.numeric(difftime(Sys.time(), begin2, units =  "secs"));
  curr_score <- score_method(curr_fit, true_betas, num_orig, num_aug, Xv, Yv)
  
  all_scores <- 
    bind_rows(all_scores, 
              bind_cols(data_seed = data_seeds[i], 
                        sim_num = i, 
                        method_name = "gim_author", 
                        run_time = curr_run_time,
                        curr_score))
  all_beta_ests <- 
    bind_rows(all_beta_ests, 
              bind_cols(data_seed = data_seeds[i], 
                        sim_num = i, 
                        method_name = "gim_author", 
                        curr_score,
                        beta_label = all_var_names,
                        true_betas = true_betas, 
                        est_betas = curr_fit$beta_hat[all_var_names], 
                        orig = unlist(map2(c(T,F), c(num_orig, num_aug), rep))))
  # Because the gim function does not return values of the lagrangian parameters
  # we estimate them here, fixing beta and sigma_h at their values. 
  if(curr_fit$converge) {
    foo <- fxn_gim_logistic_saddlepoint(Y = Yc, X = Xc,
                                        theta_tilde_with_intercept = theta_tilde_with_intercept,
                                        beta_init_with_intercept = c(curr_fit$intercept_hat,curr_fit$beta_hat),
                                        n_hist = n_hist, 
                                        Sigma_h = solve(curr_fit$Sigma_h_inv),
                                        fix_beta = TRUE,
                                        tol = convergence_tol, 
                                        max_abs_lambda = Inf)
    all_lambda_ests <- 
      bind_rows(all_lambda_ests, 
                bind_cols(data_seed =  data_seeds[i], 
                          sim_num = i, 
                          method_name = "gim_author", 
                          lambda_label = c("(Intercept)", orig_var_names),
                          est_lambda = foo$lambda_hat))
    rm(foo)
  } else {
    all_lambda_ests <- 
      bind_rows(all_lambda_ests, 
                bind_cols(data_seed =  data_seeds[i], 
                          sim_num = i, 
                          method_name = "gim_author", 
                          lambda_label = c("(Intercept)", orig_var_names),
                          est_lambda = NA_real_))
  }
  
  rm(curr_fit, curr_run_time, curr_score)
  
  # GMM (using R package GENMETA) ----
  begin2 <- Sys.time();
  curr_fit <- fxn_genmeta_logistic_author(X = Xc, 
                                          beta_internal_with_intercept = standard_mle, 
                                          theta_tilde_with_intercept = theta_tilde_with_intercept, 
                                          n_hist = n_hist, 
                                          tol = convergence_tol);
  curr_run_time <- as.numeric(difftime(Sys.time(), begin2, units =  "secs"));
  curr_score <- score_method(curr_fit, true_betas, num_orig, num_aug, Xv, Yv)
  
  all_scores <- 
    bind_rows(all_scores, 
              bind_cols(data_seed = data_seeds[i], 
                        sim_num = i, 
                        method_name = "genmeta_author", 
                        run_time = curr_run_time,
                        curr_score))
  all_beta_ests <- 
    bind_rows(all_beta_ests, 
              bind_cols(data_seed = data_seeds[i], 
                        sim_num = i, 
                        method_name = "genmeta_author", 
                        curr_score,
                        beta_label = all_var_names,
                        true_betas = true_betas, 
                        est_betas = curr_fit$beta_hat[all_var_names], 
                        orig = unlist(map2(c(T,F), c(num_orig, num_aug), rep))))
  
  # the genmeta function does not return values of the lagrangian parameters
  all_lambda_ests <- 
    bind_rows(all_lambda_ests, 
              bind_cols(data_seed =  data_seeds[i], 
                        sim_num = i, 
                        method_name = "genmeta_author", 
                        lambda_label = c("(Intercept)", orig_var_names),
                        est_lambda = NA_real_))
  
  rm(curr_fit, curr_run_time, curr_score)
  
  
  # SAB, historical value of var_theta_tilde ----
  # beta^o + beta^a ~ N(theta_tilde, eta * var_theta_tilde)
  if(!skip_bayes) {
    eigendecomp_hist_var = eigen(var_theta_tilde);
    aug_projection = create_projection(x_curr_orig = Xc_orig,
                                       x_curr_aug = Xc_aug,
                                       eigenvec_hist_var = t(eigendecomp_hist_var$vectors),
                                       imputes_list = list(c(1, 500)))[[1]];
    sab <-
      glm_sab(y = Yc, 
              x_standardized = Xc, 
              family = "binomial",
              alpha_prior_mean = theta_tilde,
              alpha_prior_cov = var_theta_tilde,
              aug_projection = aug_projection,
              phi_dist = "trunc_norm",
              phi_mean = 0.80,
              phi_sd = 0.20,
              eta_param = 2.5,
              beta_orig_scale = beta_scale_varying_phi, 
              beta_aug_scale = beta_scale_varying_phi, 
              local_dof = 1, 
              global_dof = 1, 
              slab_dof = slab_dof, 
              slab_scale = slab_scale,
              only_prior = F, 
              mc_warmup = 2e3, 
              mc_iter_after_warmup = 1e3, 
              mc_chains = 2, 
              mc_thin = 1,
              mc_stepsize = 0.1,
              mc_adapt_delta = 0.999,
              mc_max_treedepth = 14,
              eigendecomp_hist_var = eigendecomp_hist_var,
              seed = data_seeds[i]);
    curr_run_time <- as.numeric(difftime(Sys.time(), begin2, units =  "secs"));
    
    begin2 <- Sys.time();
    curr_fit <-
      list(beta_hat = apply(sab$beta, 2, median) %>% `names<-`(colnames(Xc)),
           intercept_hat = median(sab$mu) %>% `names<-`("(Intercept)"),
           converge = NA, 
           message = NA, 
           iter = NA,
           singular_hessian = NA,
           final_diff = NA,
           objective_function = c("total" = NA),
           phi = median(sab$phi), 
           eta = median(sab$eta))
    curr_score <- score_method(curr_fit, true_betas, num_orig, num_aug, Xv, Yv)
    
    all_scores <- 
      bind_rows(all_scores, 
                bind_cols(data_seed =  data_seeds[i], 
                          sim_num = i, 
                          method_name = "sab", 
                          run_time = curr_run_time,
                          curr_score))
    all_beta_ests <- 
      bind_rows(all_beta_ests, 
                bind_cols(data_seed =  data_seeds[i], 
                          sim_num = i, 
                          method_name = "sab", 
                          curr_score, 
                          beta_label = all_var_names,
                          true_betas = true_betas, 
                          est_betas = curr_fit$beta_hat[all_var_names], 
                          orig = unlist(map2(c(T,F), c(num_orig, num_aug), rep))))
    
    all_bayesian_diag <- 
      bind_rows(all_bayesian_diag, 
                bind_cols(data_seed =  data_seeds[i], 
                          sim_num = i, 
                          method_name = "sab", 
                          num_divergences = sab$num_divergences,
                          num_max_treedepth = sab$num_max_treedepth,
                          min_ebfmi = sab$min_ebfmi,
                          max_rhat = sab$max_rhat))
    
    rm(sab, curr_fit, curr_score)
    rm(aug_projection)
    
  }
  # SAB2, lm_proj, omega_sq_in_variance = TRUE ----
  # beta^o + beta^a ~ N(omega * theta_tilde, omega^2 * var_theta_tilde)
  
  if(!skip_bayes) {
    
    # Estimate the X_aug ~ X_orig multivariate regression coefficients using
    # weakly penalized ridge regression: beta ~ N(0, 25 * sigma^2)
    lambda_matrix = diag(c(0, rep(0.04, num_orig)))
    aug_projection = as.matrix(coef(lm(rbind(Xc_aug,matrix(0, nrow = num_orig + 1, ncol = num_aug)) ~ -1 + rbind(cbind(1, Xc_orig), lambda_matrix))))[-1,,drop = FALSE]
    dimnames(aug_projection) = list(orig_var_names, aug_var_names);
    rm(lambda_matrix)
    
    sab2 <-
      glm_sab2(y = Yc, 
               x_standardized = Xc, 
               family = "binomial",
               alpha_prior_mean = theta_tilde,
               alpha_prior_cov = var_theta_tilde,
               aug_projection = aug_projection,
               phi_dist = "trunc_norm",
               phi_mean = 0.80,
               phi_sd = 0.20,
               eta_param = Inf,
               omega_mean = 0, 
               omega_sd = 0.20,
               omega_sq_in_variance = TRUE,
               beta_orig_scale = beta_scale_varying_phi, 
               beta_aug_scale = beta_scale_varying_phi, 
               local_dof = 1, 
               global_dof = 1, 
               slab_dof = slab_dof, 
               slab_scale = slab_scale,
               only_prior = F, 
               mc_warmup = 2e3, 
               mc_iter_after_warmup = 1e3, 
               mc_chains = 2, 
               mc_thin = 1,
               mc_stepsize = 0.1,
               mc_adapt_delta = 0.999,
               mc_max_treedepth = 14,
               eigendecomp_hist_var = eigendecomp_hist_var,
               seed = data_seeds[i]);
    curr_run_time <- as.numeric(difftime(Sys.time(), begin2, units =  "secs"));
    
    begin2 <- Sys.time();
    curr_fit <-
      list(beta_hat = apply(sab2$beta, 2, median) %>% `names<-`(colnames(Xc)),
           intercept_hat = median(sab2$mu) %>% `names<-`("(Intercept)"),
           converge = NA, 
           message = NA, 
           iter = NA,
           singular_hessian = NA,
           final_diff = NA,
           objective_function = c("total" = NA),
           phi = median(sab2$phi), 
           eta = median(sab2$eta),
           omega = median(sab2$omega))
    curr_score <- score_method(curr_fit, true_betas, num_orig, num_aug, Xv, Yv)
    
    all_scores <- 
      bind_rows(all_scores, 
                bind_cols(data_seed =  data_seeds[i], 
                          sim_num = i, 
                          method_name = "sab2", 
                          run_time = curr_run_time,
                          curr_score))
    all_beta_ests <- 
      bind_rows(all_beta_ests, 
                bind_cols(data_seed = data_seeds[i], 
                          sim_num = i, 
                          method_name = "sab2",
                          curr_score,
                          beta_label = all_var_names,
                          true_betas = true_betas, 
                          est_betas = curr_fit$beta_hat[all_var_names], 
                          orig = unlist(map2(c(T,F), c(num_orig, num_aug), rep))))
    
    all_bayesian_diag <- 
      bind_rows(all_bayesian_diag, 
                bind_cols(data_seed =  data_seeds[i], 
                          sim_num = i, 
                          method_name = "sab2", 
                          num_divergences = sab2$num_divergences,
                          num_max_treedepth = sab2$num_max_treedepth,
                          min_ebfmi = sab2$min_ebfmi,
                          max_rhat = sab2$max_rhat))
    
    rm(sab2, curr_fit, curr_score)
  }
  
  # SAB2, lm_proj, omega_sq_in_variance = FALSE ----
  # beta^o + beta^a ~ N(omega * theta_tilde, eta * var_theta_tilde)
  if(!skip_bayes) {
    sab2 <-
      glm_sab2(y = Yc, 
               x_standardized = Xc, 
               family = "binomial",
               alpha_prior_mean = theta_tilde,
               alpha_prior_cov = var_theta_tilde,
               aug_projection = aug_projection,
               phi_dist = "trunc_norm",
               phi_mean = 0.80,
               phi_sd = 0.20,
               eta_param = 2.5, # 1st distinguishing feature from above version of SAB2
               omega_mean = 0, 
               omega_sd = 0.20,
               omega_sq_in_variance = FALSE, # 2st distinguishing feature from above version of SAB2
               beta_orig_scale = beta_scale_varying_phi, 
               beta_aug_scale = beta_scale_varying_phi, 
               local_dof = 1, 
               global_dof = 1, 
               slab_dof = slab_dof, 
               slab_scale = slab_scale,
               only_prior = F, 
               mc_warmup = 2e3, 
               mc_iter_after_warmup = 1e3, 
               mc_chains = 2, 
               mc_thin = 1,
               mc_stepsize = 0.1,
               mc_adapt_delta = 0.999,
               mc_max_treedepth = 14,
               eigendecomp_hist_var = eigendecomp_hist_var,
               seed = data_seeds[i]);
    curr_run_time <- as.numeric(difftime(Sys.time(), begin2, units =  "secs"));
    
    begin2 <- Sys.time();
    curr_fit <-
      list(beta_hat = apply(sab2$beta, 2, median) %>% `names<-`(colnames(Xc)),
           intercept_hat = median(sab2$mu) %>% `names<-`("(Intercept)"),
           converge = NA, 
           message = NA, 
           iter = NA,
           singular_hessian = NA,
           final_diff = NA,
           objective_function = c("total" = NA),
           phi = median(sab2$phi), 
           eta = median(sab2$eta),
           omega = median(sab2$omega))
    curr_score <- score_method(curr_fit, true_betas, num_orig, num_aug, Xv, Yv)
    
    all_scores <- 
      bind_rows(all_scores, 
                bind_cols(data_seed = data_seeds[i], 
                          sim_num = i,
                          method_name = "sab2_alt", 
                          run_time = curr_run_time,
                          curr_score))
    all_beta_ests <- 
      bind_rows(all_beta_ests, 
                bind_cols(data_seed = data_seeds[i], 
                          sim_num = i, 
                          method_name = "sab2_alt",
                          curr_score,
                          beta_label = all_var_names,
                          true_betas = true_betas, 
                          est_betas = curr_fit$beta_hat[all_var_names], 
                          orig = unlist(map2(c(T,F), c(num_orig, num_aug), rep))))
    
    all_bayesian_diag <- 
      bind_rows(all_bayesian_diag, 
                bind_cols(data_seed =  data_seeds[i], 
                          sim_num = i, 
                          method_name = "sab2_alt", 
                          num_divergences = sab2$num_divergences,
                          num_max_treedepth = sab2$num_max_treedepth,
                          min_ebfmi = sab2$min_ebfmi,
                          max_rhat = sab2$max_rhat))
    
    rm(sab2, curr_fit, curr_score)
    rm(eigendecomp_hist_var)
  }
  rm(var_theta_tilde);
  
  # SAB2, lm_proj, omega_sq_in_variance = TRUE, internal estimate of var_theta ----
  # beta^o + beta^a ~ N(omega * theta_tilde, omega^2 * var_theta_tilde)
  if(!skip_bayes) {
    
    var_theta_tilde_alt <- (n_curr / n_hist) * vcov(logistf(Yc ~ Xc_orig, family="binomial", plconf = 1))[-1, -1, drop = FALSE]
    
    sab2 <-
      glm_sab2(y = Yc, 
               x_standardized = Xc, 
               family = "binomial",
               alpha_prior_mean = theta_tilde,
               alpha_prior_cov = var_theta_tilde_alt,
               aug_projection = aug_projection,
               phi_dist = "trunc_norm",
               phi_mean = 0.80,
               phi_sd = 0.20,
               eta_param = Inf,
               omega_mean = 0, 
               omega_sd = 0.20,
               omega_sq_in_variance = TRUE,
               beta_orig_scale = beta_scale_varying_phi, 
               beta_aug_scale = beta_scale_varying_phi, 
               local_dof = 1, 
               global_dof = 1, 
               slab_dof = slab_dof, 
               slab_scale = slab_scale,
               only_prior = F, 
               mc_warmup = 2e3, 
               mc_iter_after_warmup = 1e3, 
               mc_chains = 2, 
               mc_thin = 1,
               mc_stepsize = 0.1,
               mc_adapt_delta = 0.999,
               mc_max_treedepth = 14,
               eigendecomp_hist_var = NULL,
               seed = data_seeds[i]);
    curr_run_time <- as.numeric(difftime(Sys.time(), begin2, units =  "secs"));
    
    begin2 <- Sys.time();
    curr_fit <-
      list(beta_hat = apply(sab2$beta, 2, median) %>% `names<-`(colnames(Xc)),
           intercept_hat = median(sab2$mu) %>% `names<-`("(Intercept)"),
           converge = NA, 
           message = NA,   
           iter = NA,
           singular_hessian = NA,
           final_diff = NA,
           objective_function = c("total" = NA),
           phi = median(sab2$phi), 
           eta = median(sab2$eta),
           omega = median(sab2$omega))
    curr_score <- score_method(curr_fit, true_betas, num_orig, num_aug, Xv, Yv)
    
    all_scores <- 
      bind_rows(all_scores, 
                bind_cols(data_seed =  data_seeds[i],
                          sim_num = i,
                          method_name = "sab2_flexible", 
                          run_time = curr_run_time,
                          curr_score))
    all_beta_ests <- 
      bind_rows(all_beta_ests, 
                bind_cols(data_seed =  data_seeds[i], 
                          sim_num = i, 
                          method_name = "sab2_flexible",
                          curr_score,
                          beta_label = all_var_names,
                          true_betas = true_betas, 
                          est_betas = curr_fit$beta_hat[all_var_names], 
                          orig = unlist(map2(c(T,F), c(num_orig, num_aug), rep))))
    
    all_bayesian_diag <- 
      bind_rows(all_bayesian_diag, 
                bind_cols(data_seed =  data_seeds[i], 
                          sim_num = i, 
                          method_name = "sab2_flexible", 
                          num_divergences = sab2$num_divergences,
                          num_max_treedepth = sab2$num_max_treedepth,
                          min_ebfmi = sab2$min_ebfmi,
                          max_rhat = sab2$max_rhat))
    
    rm(sab2, curr_fit, curr_score)
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
  
  rm(standard_mle)
  
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
  
  all_beta_ests <- 
    all_beta_ests %>%
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
  
  all_lambda_ests <- 
    all_lambda_ests %>%
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
           method_name,
           everything())
  
  all_bayesian_diag <- 
    all_bayesian_diag %>%
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
           method_name,
           everything())
  
  if(!my_computer) {
    # Save scores and estimates and running times as csv
    write_csv(all_scores,
              file = paste0("out/job",array_id_with_offset, "_scores.csv"),
              append = FALSE);
    
    write_csv(all_beta_ests,
              file = paste0("out/job",array_id_with_offset, "_betas.csv"),
              append = FALSE);
    
    write_csv(all_lambda_ests,
              file = paste0("out/job",array_id_with_offset, "_lambdas.csv"),
              append = FALSE);
    
    write_csv(all_bayesian_diag,
              file = paste0("out/job",array_id_with_offset, "_bayesian_diag.csv"),
              append = FALSE);
    
    write_csv(array_id_stats,
              file = paste0("out/job",array_id_with_offset, "_array_stats.csv"),
              append = FALSE);
    
  }
}

if(my_computer) {
  all_beta_ests %>% 
    filter(beta_label == "p1") %>%
    select(sim_num, method_name, est_betas) %>% 
    pivot_wider(id_cols = sim_num, names_from = method_name, values_from = est_betas)
  all_scores %>% 
    group_by(method_name) %>% 
    summarize(mean(squared_error), mean(squared_error_orig), mean(squared_error_aug))
}

