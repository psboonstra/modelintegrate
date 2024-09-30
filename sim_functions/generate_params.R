library(tidyverse)
# the effect sizes are select using the Rscript OR_from_AUC.R to find the
# combination that ensures a specific AUC

scenario_list <- 
  list(
    #auc_to_betas(0.85, relative_scales = c(1, 1, 2), cor = 0.25, which_binary = c(F, F, T), seed = 1)$betas;
    scenario1 =
      list(true_bivariate_cor = 0.25, 
           orig = c(0.594, 0.594),
           aug = 1.189, 
           which_binary = c(F, F, T)),
    #auc_to_betas(0.75, relative_scales = c(4, 1, -4, -4, -1, -1, 1, 1), cor = 0.60, which_binary = c(F, F, T, T, T, T, T, T), seed = 1)$betas;
    scenario2 = 
      list(true_bivariate_cor = 0.60, 
           orig = c(0.736, 0.184), 
           aug = c(-0.736, -0.736, -0.184, -0.184, 0.184, 0.184),
           which_binary = c(F, F, T, T, T, T, T, T)),
    #auc_to_betas(0.7, relative_scales = rep(1, 30), cor = 0.40, which_binary = rep(T, 30), seed = 1)$betas;
    scenario3 = 
      list(true_bivariate_cor = 0.40, 
           orig = c(0.047, 0.047, 0.047, 0.047, 0.047,
                    0.047, 0.047, 0.047, 0.047, 0.047, 
                    0.047, 0.047, 0.047, 0.047, 0.047, 
                    0.047, 0.047, 0.047, 0.047, 0.047),
           aug = c(0.047, 0.047, 0.047, 0.047, 0.047, 
                   0.047, 0.047, 0.047, 0.047, 0.047),
           which_binary = rep(T, 30)),
    #auc_to_betas(0.7, relative_scales = rep(c(20, -1, -1, -1, -1, -1, 10, 10), each = 5),cor = 0.025, which_binary = c(rep(F, 30), rep(T, 10)), seed = 1)$betas;
    scenario4 = 
      list(true_bivariate_cor = 0.025, 
           orig = c(0.262, 0.262, 0.262, 0.262, 0.262,
                    -0.013, -0.013, -0.013, -0.013, -0.013, 
                    -0.013, -0.013, -0.013, -0.013, -0.013, 
                    -0.013, -0.013, -0.013, -0.013, -0.013, 
                    -0.013, -0.013, -0.013, -0.013, -0.013, 
                    -0.013, -0.013, -0.013, -0.013, -0.013), 
           aug = c(0.131, 0.131, 0.131, 0.131, 0.131, 
                   0.131, 0.131, 0.131, 0.131, 0.131),
           which_binary = c(rep(F, 30), rep(T, 10)))
  )

scenario_num_orig <-
  map(scenario_list, ~.x[["orig"]]) %>%
  map_int(length)

true_beta0 <- -1