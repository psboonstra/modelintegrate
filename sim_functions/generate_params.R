library(tidyverse)
# the effect sizes are select using the Rscript OR_from_AUC.R to find the combination that ensures a specific AUC
# source("OR_from_AUC.R")


scenario_list <- 
  list(
    #auc_to_betas(0.85, relative_scales = c(1, 1, 2), cor = 0.25)$betas;
    scenario1 =
      list(true_bivariate_cor = 0.25, 
           orig = c(0.637, 0.637),
           aug = 1.274),
    #auc_to_betas(0.75, relative_scales = c(25, 10, 1, 1, -15, 15), 0.10)$betas;
    scenario2 = 
      list(true_bivariate_cor = 0.10, 
           orig = c(0.771, 0.309, 0.031, 0.031), 
           aug = c(-0.463, 0.463)),
    #auc_to_betas(0.7, relative_scales = rep(1, 30), 0.05)$betas;
    scenario3 = 
      list(true_bivariate_cor = 0.05, 
           orig = c(0.092, 0.092, 0.092, 0.092, 0.092,
                    0.092, 0.092, 0.092, 0.092, 0.092, 
                    0.092, 0.092, 0.092, 0.092, 0.092, 
                    0.092, 0.092, 0.092, 0.092, 0.092),
           aug = c(0.092, 0.092, 0.092, 0.092, 0.092, 
                   0.092, 0.092, 0.092, 0.092, 0.092)),
    #auc_to_betas(0.7, relative_scales = rep(c(20, -1, -1, -1, -1, -1, 10, 10), each = 5), 0.025)$betas;
    scenario4 = 
      list(true_bivariate_cor = 0.025, 
           orig = c(0.261, 0.261, 0.261, 0.261, 0.261,
                    -0.013, -0.013, -0.013, -0.013, -0.013, 
                    -0.013, -0.013, -0.013, -0.013, -0.013, 
                    -0.013, -0.013, -0.013, -0.013, -0.013, 
                    -0.013, -0.013, -0.013, -0.013, -0.013, 
                    -0.013, -0.013, -0.013, -0.013, -0.013), 
           aug = c(0.130, 0.130, 0.130, 0.130, 0.130, 
                   0.130, 0.130, 0.130, 0.130, 0.130))
  )

scenario_num_orig <-
  map(scenario_list, ~.x[["orig"]]) %>%
  map_int(length)

true_beta0 <- -1