## ----message = F, warning = F----
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
library(logistf);
library(knitr);
library(kableExtra);


## -----------------------------
sims_per_scenario = 200;
which_batch = 1;
total_num_arrays_by_batch = c(999, 999, 999)
total_num_arrays = sum(total_num_arrays_by_batch)


## -----------------------------
source("sim_functions/OR_from_AUC.R")
source("sim_functions/generate_params.R")
source("sim_functions/scenario_setup.R")
source("sim_functions/score_method.R")
source("aux_functions/log1plex.R")
source("aux_functions/get_num_eff.R")


## -----------------------------
col_types_scores <- cols(
  array_id = col_integer(), 
  sim_num = col_integer(), 
  which_batch = col_integer(), 
  data_seed = col_integer(), 
  n_hist = col_integer(), 
  n_curr = col_integer(), 
  true_bivariate_cor = col_double(), 
  different_external_model = col_logical(),
  different_covariate_dist = col_logical(),
  scenario_name = col_character(),
  method_name = col_character(),
  squared_error = col_double(), 
  squared_error_orig = col_double(), 
  squared_error_aug = col_double(), 
  brier_score = col_double(),
  auc = col_double(),
  converge = col_logical(), 
  message = col_character(), 
  iter = col_integer(), 
  singular = col_logical(), 
  final_diff = col_double(),
  final_objective_function = col_double()
);

col_types_betas <- cols(
  array_id = col_integer(), 
  sim_num = col_integer(), 
  which_batch = col_integer(), 
  data_seed = col_integer(), 
  n_hist = col_integer(), 
  n_curr = col_integer(), 
  true_bivariate_cor = col_double(), 
  different_external_model = col_logical(),
  different_covariate_dist = col_logical(),
  scenario_name = col_character(),
  method_name = col_character(),
  true_betas = col_double(), 
  orig = col_logical(), 
  beta_label = col_character(),
  est_betas = col_double(),
  squared_error = col_double(), 
  squared_error_orig = col_double(), 
  squared_error_aug = col_double(), 
  brier_score = col_double(),
  auc = col_double(),
  converge = col_logical(), 
  message = col_character(), 
  iter = col_integer(), 
  singular = col_logical(), 
  final_diff = col_double(),
  final_objective_function = col_double()
);

col_types_lambdas <- cols(
  array_id = col_integer(), 
  sim_num = col_integer(), 
  which_batch = col_integer(), 
  data_seed = col_integer(), 
  n_hist = col_integer(), 
  n_curr = col_integer(), 
  true_bivariate_cor = col_double(), 
  different_external_model = col_logical(),
  different_covariate_dist = col_logical(),
  scenario_name = col_character(),
  method_name = col_character(),
  lambda_label = col_character(),
  est_lambda = col_double(),
);

col_types_array_stats <- cols(
  n_hist = col_double(),
  n_curr = col_double(),
  different_external_model = col_logical(),
  different_covariate_dist = col_logical(),
  scenario_name = col_character(),
  scenario = col_double(),
  num_orig = col_double(),
  num_aug = col_double(),
  old_scenario_complexity_weight = col_double(),
  scenario_complexity_weight = col_double(),
  num_arrays = col_double(),
  start_array_id = col_double(),
  end_array_id = col_double(),
  array_id = col_double(),
  n_sim_this_array_id = col_double(),
  total_runtime_secs = col_double()
);

col_types_bayesian_diag <- cols(
  array_id = col_integer(), 
  sim_num = col_integer(), 
  which_batch = col_integer(), 
  data_seed = col_integer(), 
  n_hist = col_integer(), 
  n_curr = col_integer(), 
  true_bivariate_cor = col_double(), 
  different_external_model = col_logical(),
  different_covariate_dist = col_logical(),
  scenario_name = col_character(),
  method_name = col_character(),
  num_divergences = col_integer(),
  num_max_treedepth = col_integer(), 
  min_ebfmi = col_double(), 
  max_rhat = col_double()
);


## -----------------------------
if(!exists("read_all_sims")) {read_all_sims = TRUE;}
if(read_all_sims) {
  which_sims = 1:total_num_arrays
} else {
  which_sims = sample(total_num_arrays, 50);
}
not_found = NULL


for(i in which_sims) {
  foo <- try(read_csv(paste0("out/job",i,"_scores.csv"), col_types = col_types_scores));
  if(!"try-error" %in% class(foo)) {
    assign(paste0("job",i,"_scores"), foo)
  } else {
    cat("sim ", i, ", not found\n");
    not_found = c(not_found,i)
  }
  foo <- try(read_csv(paste0("out/job",i,"_betas.csv"), col_types = col_types_betas));
  if(!"try-error" %in% class(foo)) {
    assign(paste0("job",i,"_betas"), foo)
  } 
  foo <- try(read_csv(paste0("out/job",i,"_lambdas.csv"), col_types = col_types_lambdas));
  if(!"try-error" %in% class(foo)) {
    assign(paste0("job",i,"_lambdas"), foo)
  } 
  foo <- try(read_csv(paste0("out/job",i,"_bayesian_diag.csv"), col_types = col_types_bayesian_diag));
  if(!"try-error" %in% class(foo)) {
    assign(paste0("job",i,"_bayesian_diag"), foo)
  } 
  #if(foo %>% filter(data_seed==2080255079, method_name == "sab2_flexible", num_divergences==316) %>% nrow()) {stop();}
  foo <- try(read_csv(paste0("out/job",i,"_array_stats.csv"), col_types = col_types_array_stats));
  if(!"try-error" %in% class(foo)) {
    assign(paste0("job",i,"_array_stats"), foo)
  } 
  
  if(i%%500 == 0) {cat(i, "\n")}
}

rm(col_types_betas, col_types_scores, col_types_lambdas, col_types_bayesian_diag, col_types_array_stats)

raw_all_scores <-
  bind_rows(map(str_subset(ls(), pattern = "job\\d+_scores"), get)) %>%
  arrange(array_id, sim_num)
raw_all_betas <-
  bind_rows(map(str_subset(ls(), pattern = "job\\d+_betas"), get)) %>%
  arrange(array_id, sim_num)
raw_all_lambdas <-
  bind_rows(map(str_subset(ls(), pattern = "job\\d+_lambdas"), get)) %>%
  arrange(array_id, sim_num)
raw_all_bayesian_diag <-
  bind_rows(map(str_subset(ls(), pattern = "job\\d+_bayesian_diag"), get)) %>%
  arrange(array_id, sim_num)
corrected_array_ids <- 
  (str_subset(ls(), pattern = "job\\d+_array_stats") %>%
     str_extract("\\d+") %>% 
     as.numeric())
raw_all_array_stats <-
  bind_rows(map(str_subset(ls(), pattern = "job\\d+_array_stats"), get)) %>%
  select(-array_id) %>%
  bind_cols(tibble(array_id = corrected_array_ids)) %>%
  arrange(array_id)
unique(raw_all_scores$method_name)


## -----------------------------
raw_all_array_stats %>% 
  group_by(n_hist, n_curr, different_external_model, different_covariate_dist,scenario_name, scenario, start_array_id, end_array_id) %>% 
  summarize(n_sim = sum(n_sim_this_array_id), total_runtime_secs = sum(total_runtime_secs),
            .groups = "drop") %>% 
  left_join(all_scenarios %>% select(scenario, num_orig, num_aug, scenario_complexity_weight, old_scenario_complexity_weight)) %>%
  ggplot() + 
  geom_point(aes(x = total_runtime_secs, 
                 y = scenario_complexity_weight, 
                 color = scenario_name, 
                 shape = factor(different_covariate_dist))) 

raw_all_array_stats %>% 
  group_by(n_hist, n_curr, different_external_model, different_covariate_dist,scenario_name, scenario, start_array_id, end_array_id) %>% 
  summarize(n_sim = sum(n_sim_this_array_id), total_runtime_secs = sum(total_runtime_secs),
            .groups = "drop") %>% 
  left_join(all_scenarios %>% select(scenario, num_orig, num_aug, scenario_complexity_weight)) %>%
  lm(formula = log(total_runtime_secs) ~ factor(scenario_name) + log(n_hist) * log(n_curr) ) %>%
  predict()


## -----------------------------

primary_method_names = c(
  "cml_saddlepoint", 
  #"gim_sandwich", 
  "gim_author",
  #"genmeta_author",
  "ratios", 
  "sab", # eta scales variance
  "sab2"#, # omega scales mean, omega^2 scales variance, linear projection matrix
)

secondary_method_names = c(
  "cml_newtonraphson", 
  "sab2_flexible" # like sab2, but uses internal covariance matrix 
)

name_mapping = 
  c("truth" = "Truth",
    "glm_vanilla" = "GLM(MLE)",
    "glm_bayes" = "GLM(Curr)",
    "cml_saddlepoint" = "CML(SP)",
    "cml_newtonraphson" = "CML(NR)",
    "cml_newtonraphson_step" = "CML(NRs)",
    "gim_sandwich" = "GIM(Swch)",
    "gim_author" = "GIM",
    "genmeta_author" = "GMM(Aut)",
    "ratios" = "Ratios",
    "sab" = "SAB",
    "sab2" = "SAB2",
    "sab2_alt" = "SAB2(Alt)",
    "sab2_flexible" = "SAB2(Ext)")



## -----------------------------
raw_all_scores %>% 
  filter(method_name == first(method_name)) %>% 
  group_by(n_hist, n_curr, true_bivariate_cor, method_name, scenario_name, 
           different_external_model, different_covariate_dist) %>% 
  summarize(n1 = n_distinct(data_seed), 
            n2 = n(), 
            .groups = "drop") %>% 
  mutate(diff = n1 - n2) %>% 
  pull(diff) %>%
  "=="(0) %>%
  all()


## -----------------------------
all_scores <- 
  full_join(
    raw_all_scores %>%
      mutate(converge =
               case_when(
                 str_detect(method_name, "bayes") ~ TRUE, 
                 str_detect(method_name, "sab") ~ TRUE, 
                 str_detect(method_name, "truth") ~ TRUE, 
                 str_detect(method_name, "genmeta_author") & is.na(converge) ~ TRUE, 
                 TRUE ~ converge
               )) %>%
      mutate(squared_error = ifelse(converge, squared_error, NA_real_),
             squared_error_orig = ifelse(converge, squared_error_orig, NA_real_),
             squared_error_aug = ifelse(converge, squared_error_aug, NA_real_)),
    raw_all_scores %>%
      filter(method_name == "glm_bayes") %>% 
      rename_with(~ paste0(.x, "_glm_bayes"), c(squared_error, squared_error_orig, squared_error_aug, converge)) %>%
      dplyr::select(array_id:scenario_name, squared_error_glm_bayes: squared_error_aug_glm_bayes, converge_glm_bayes)) %>%
  mutate(rel_squared_error = squared_error / squared_error_glm_bayes, 
         log10_rel_squared_error = log10(rel_squared_error),
         log10_squared_error = log10(squared_error),
         rel_squared_error_orig = squared_error_orig / squared_error_orig_glm_bayes, 
         log10_rel_squared_error_orig = log10(rel_squared_error_orig),
         log10_squared_error_orig = log10(squared_error_orig),
         is_primary = method_name %in% primary_method_names, 
         is_secondary = method_name %in% secondary_method_names) %>%
  # whether to group by n_hist in this next line is important:
  # - if you **do** group by n_hist, then each value of n_hist will potentially 
  # have a different set of data_seeds for which all methods converge. Do this
  # if you want to drop fewer data_seeds
  # - if you **don't** group by n_hist, then you will end up with a common set of
  # data_seeds across all values of n_hist. Do this if you have one or more methods
  # the performance of which does not depend upon n_hist. 
  group_by(n_hist, n_curr, true_bivariate_cor, scenario_name, 
           different_external_model, different_covariate_dist, data_seed) %>%
  #group_by(n_curr, true_bivariate_cor, scenario_name, 
  #         different_external_model, different_covariate_dist, data_seed) %>%
  mutate(all_primary_converge = all(converge | !is_primary)) %>%
  mutate(all_primary_secondary_converge = all(converge | (!is_primary & !is_secondary))) %>%
  ungroup() %>%
  mutate(n_curr_fancy = glue("n=={n_curr}")) %>%
  mutate(n_hist = factor(n_hist))

all_betas <- 
  raw_all_betas %>%
  mutate(converge =
           case_when(
             str_detect(method_name, "bayes") ~ TRUE, 
             str_detect(method_name, "sab") ~ TRUE, 
             str_detect(method_name, "truth") ~ TRUE, 
             str_detect(method_name, "genmeta_author") & is.na(converge) ~ TRUE, 
             TRUE ~ converge
           ),
         est_betas = ifelse(converge, est_betas, NA_real_),
         is_primary = method_name %in% primary_method_names,
         is_secondary = method_name %in% secondary_method_names) %>%
  # whether to group by n_hist in this next line is important:
  # - if you **do** group by n_hist, then each value of n_hist will potentially 
  # have a different set of data_seeds for which all methods converge. Do this
  # if you want to drop fewer data_seeds
  # - if you **don't** group by n_hist, then you will end up with a common set of
  # data_seeds across all values of n_hist. Do this if you have one or more methods
  # the performance of which does not depend upon n_hist. 
  group_by(n_hist, n_curr, true_bivariate_cor, scenario_name, 
           different_external_model, different_covariate_dist, data_seed) %>%
  # group_by(n_curr, true_bivariate_cor, scenario_name, 
  #         different_external_model, different_covariate_dist, data_seed) %>%
  mutate(all_primary_converge = all(converge | !is_primary)) %>%
  mutate(all_primary_secondary_converge = all(converge | (!is_primary & !is_secondary))) %>%
  ungroup() %>%
  mutate(n_curr_fancy = glue("n=={n_curr}")) %>%
  mutate(n_hist = factor(n_hist))


## -----------------------------
rimse_figure_data <- 
  all_scores %>% 
  filter(all_primary_converge, 
         is_primary) %>%
  group_by(scenario_name, true_bivariate_cor, 
           n_curr, n_hist, n_curr_fancy, 
           method_name,
           different_external_model, different_covariate_dist) %>%
  mutate(median_log10_rel_squared_error_orig = median(log10_rel_squared_error_orig), 
         num_datasets = n()) %>%
  ungroup(method_name) %>%
  mutate(rank_log10_rel_squared_error_orig = rank(median_log10_rel_squared_error_orig, ties = "max") / num_datasets) %>%
  ungroup() %>%
  mutate(num_datasets = glue("{num_datasets} / {max(num_datasets)}") %>% as.character(),
         method_name = 
           as.character(name_mapping[method_name]) %>% 
           factor(levels = name_mapping[primary_method_names]),
         base_case = !different_external_model & !different_covariate_dist,
         setting = 
           case_when(!different_external_model & !different_covariate_dist ~ "No differences",
                     different_external_model & !different_covariate_dist ~ "Covariate model",
                     !different_external_model & different_covariate_dist ~ "Outcome model",
                     different_external_model & different_covariate_dist ~ "Covariate and outcome model")) %>%
  arrange(scenario_name, true_bivariate_cor, 
          n_curr, n_hist, 
          method_name, 
          different_external_model, different_covariate_dist) %>% 
  mutate(scenario_name_fancy = 
           case_when(
             scenario_name == "scenario1" ~ "list(beta^o==group('{',list(0.594,0.594),'}'),beta^a==1.189)",
             scenario_name == "scenario2" ~ "list(beta^o==group('{',list(0.736, 0.184),'}'),beta^a==group('{',list(-0.736, -0.736, -0.184, -0.184, 0.184, 0.184),'}'))",
             scenario_name == "scenario3" ~ "list(group('{',list(p,q),'}')==group('{',list(20,10),'}'), uniform~signal)",
             scenario_name == "scenario4" ~ "list(group('{',list(p,q),'}')==group('{',list(30,10),'}'), nonuniform~signal)"
           ) %>% fct_inorder(),
         scenario_name_fancy_simple = 
           case_when(
             scenario_name == "scenario1" ~ glue("Scen~1~':'~group('{{',list(p,q),'}}')==group('{{',list(2,1),'}}')~';'~{n_curr_fancy}"),
             scenario_name == "scenario2" ~ glue("Scen~2~':'~group('{{',list(p,q),'}}')==group('{{',list(2,6),'}}')~';'~{n_curr_fancy}"),
             scenario_name == "scenario3" ~ glue("Scen~3~':'~group('{{',list(p,q),'}}')==group('{{',list(20,10),'}}')~';'~{n_curr_fancy}"),
             scenario_name == "scenario4" ~ glue("Scen~4~':'~group('{{',list(p,q),'}}')==group('{{',list(30,10),'}}')~';'~{n_curr_fancy}")
           ) %>% fct_inorder()) %>%
  group_by(n_hist, n_curr, n_curr_fancy, true_bivariate_cor, different_external_model, different_covariate_dist, 
           scenario_name, scenario_name_fancy_simple, base_case, 
           method_name, 
           num_datasets, rank_log10_rel_squared_error_orig) %>%
  summarize(y0 = quantile(log10_rel_squared_error_orig, 0.25),
            y25 = quantile(log10_rel_squared_error_orig, 0.25), 
            y50 = quantile(log10_rel_squared_error_orig, 0.5), 
            y75 = quantile(log10_rel_squared_error_orig, 0.75), 
            y100 = quantile(log10_rel_squared_error_orig, 0.75), 
            .groups = "drop") 


filter(rimse_figure_data, base_case) %>%
  ggplot(mapping = aes(x = n_hist)) + 
  facet_wrap(facets = vars(scenario_name_fancy_simple), 
             nrow = 4,
             labeller = label_parsed, scales = "free_y") + 
  scale_fill_manual(name = "Method", 
                    values = RColorBrewer::brewer.pal(length(primary_method_names), "Dark2")) + 
  scale_y_continuous(name = "logRIMSE", expand = expansion(mult = c(0.15, 0.15))) + 
  scale_x_discrete(name = expression(n[h]), expand = expansion(mult = 0.15)) + 
  theme(
    text = element_text(size = 18), 
    strip.text = element_text(size = 14),
    legend.key.width = unit(1., "cm"),
    legend.position = "top") + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_boxplot(mapping = aes(ymin = y0, 
                             lower = y25, 
                             middle = y50,
                             upper = y75, 
                             ymax = y100, 
                             fill = method_name),
               stat = "identity") + 
  geom_text(mapping = aes(label = rank_log10_rel_squared_error_orig,
                          group = method_name, 
                          y = -Inf),
            vjust = -0.2,
            position = position_dodge(width = 0.9)) +
  geom_text(data = rimse_figure_data %>% 
              filter(base_case) %>% 
              group_by(scenario_name, true_bivariate_cor, 
                       n_curr, n_hist, n_curr_fancy, 
                       different_external_model, different_covariate_dist) %>%
              slice(1),
            mapping = aes(label = num_datasets,
                          y = Inf),
            vjust = 1.2) +
  guides(fill = guide_legend(nrow = 1))
ggsave(filename = "figures/figure1.png", height = 12, width = 10, dpi = 500)
ggsave(filename = "figures/figure1.eps", height = 12, width = 10)


rimse_figure_data %>% 
  filter(base_case | (!different_external_model & different_covariate_dist)) %>%
  ggplot(mapping = aes(x = n_hist)) + 
  facet_wrap(facets = vars(scenario_name_fancy_simple), 
             nrow = 4,
             labeller = label_parsed, scales = "free_y") + 
  scale_fill_manual(name = "Method", 
                    values = RColorBrewer::brewer.pal(length(primary_method_names), "Dark2")) + 
  scale_y_continuous(name = "logRIMSE", expand = expansion(mult = c(0.15, 0.025))) + 
  scale_x_discrete(name = expression(n[h]), expand = expansion(mult = 0.15)) +
  scale_alpha_manual(name = "Differences in\ngenerating models", 
                     values = c(1, 0.),
                     breaks = c(FALSE, TRUE), 
                     labels = c(expression("Different"~'['~list(X^o, X^a)~']'), "No differences")) +
  scale_shape_manual(name = "Differences in\ngenerating models", 
                     values = c(8, NA),
                     breaks = c(FALSE, TRUE), 
                     labels = c(expression("Different"~'['~list(X^o, X^a)~']'), "No differences")) +
  theme(text = element_text(size = 18), 
        strip.text = element_text(size = 14),
        legend.box = "vertical",
        legend.key.width = unit(1., "cm"),
        legend.position = "top") + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_boxplot(mapping = aes(ymin = y0, 
                             lower = y25, 
                             middle = y50,
                             upper = y75, 
                             ymax = y100, 
                             fill = method_name,
                             alpha = base_case, 
                             group = interaction(base_case, method_name, n_hist)),
               stat = "identity") + 
  geom_point(mapping = aes(y = y50, group = interaction(base_case, method_name, n_hist), shape = base_case),
             position = position_dodge(width = 0.9)) +
  geom_text(data = rimse_figure_data %>% 
              filter(!different_external_model & different_covariate_dist),
            mapping = aes(label = rank_log10_rel_squared_error_orig,
                          group = method_name, 
                          y = -Inf),
            vjust = -0.2,
            position = position_dodge(width = 0.9)) +
  guides(fill = guide_legend(nrow = 1, order = 0),
         shape = guide_legend(nrow = 1, order = 1),
         alpha = guide_legend(nrow = 1, override.aes = list(fill = "darkgrey"), order = 1))
ggsave(filename = "figures/figure2.png", height = 12, width = 10, dpi = 500)
ggsave(filename = "figures/figure2.eps", height = 12, width = 10)


rimse_figure_data %>% 
  filter(base_case | (different_external_model & !different_covariate_dist)) %>%
  ggplot(mapping = aes(x = n_hist)) + 
  facet_wrap(facets = vars(scenario_name_fancy_simple), 
             nrow = 4,
             labeller = label_parsed, scales = "free_y") + 
  scale_fill_manual(name = "Method", 
                    values = RColorBrewer::brewer.pal(length(primary_method_names), "Dark2")) + 
  scale_y_continuous(name = "logRIMSE", expand = expansion(mult = c(0.15, 0.025))) + 
  scale_x_discrete(name = expression(n[h]), expand = expansion(mult = 0.15)) +
  scale_alpha_manual(name = "Differences in\ngenerating models", 
                     values = c(1, 0.),
                     breaks = c(FALSE, TRUE), 
                     labels = c(expression("Different"~'['~Y~"|"~list(X^o, X^a)~']'), "No differences")) +  
  scale_shape_manual(name = "Differences in\ngenerating models", 
                     values = c(8, NA),
                     breaks = c(FALSE, TRUE), 
                     labels = c(expression("Different"~'['~Y~"|"~list(X^o, X^a)~']'), "No differences")) +  
  theme(text = element_text(size = 18), 
        strip.text = element_text(size = 14),
        legend.box = "vertical",
        legend.key.width = unit(1., "cm"),
        legend.position = "top") + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_boxplot(mapping = aes(ymin = y0, 
                             lower = y25, 
                             middle = y50,
                             upper = y75, 
                             ymax = y100, 
                             fill = method_name,
                             alpha = base_case, 
                             group = interaction(base_case, method_name, n_hist)),
               stat = "identity") + 
  geom_point(mapping = aes(y = y50, group = interaction(base_case, method_name, n_hist), shape = base_case),
             position = position_dodge(width = 0.9)) +
  geom_text(data = rimse_figure_data %>% 
              filter(!different_external_model & different_covariate_dist),
            mapping = aes(label = rank_log10_rel_squared_error_orig,
                          group = method_name, 
                          y = -Inf),
            vjust = -0.2,
            position = position_dodge(width = 0.9)) +
  guides(fill = guide_legend(nrow = 1, order = 0),
         shape = guide_legend(nrow = 1, order = 1),
         alpha = guide_legend(nrow = 1, override.aes = list(fill = "darkgrey"), order = 1))
ggsave(filename = "figures/figure3.png", height = 12, width = 10, dpi = 500)
ggsave(filename = "figures/figure3.eps", height = 12, width = 10)


rimse_figure_data %>% 
  filter(base_case | (different_external_model & different_covariate_dist)) %>%
  ggplot(mapping = aes(x = n_hist)) + 
  facet_wrap(facets = vars(scenario_name_fancy_simple), 
             nrow = 4,
             labeller = label_parsed, scales = "free_y") + 
  scale_fill_manual(name = "Method", 
                    values = RColorBrewer::brewer.pal(length(primary_method_names), "Dark2")) + 
  scale_y_continuous(name = "logRIMSE", expand = expansion(mult = c(0.15, 0.025))) + 
  scale_x_discrete(name = expression(n[h]), expand = expansion(mult = 0.15)) +
  scale_alpha_manual(name = "Differences in\ngenerating models", 
                     values = c(1, 0.),
                     breaks = c(FALSE, TRUE), 
                     labels = c(expression("Different"~'['~Y~"|"~list(X^o, X^a)~']'~" and "~'['~list(X^o, X^a)~']'), "No differences")) +
  scale_shape_manual(name = "Differences in\ngenerating models", 
                     values = c(8, NA),
                     breaks = c(FALSE, TRUE), 
                     labels = c(expression("Different"~'['~Y~"|"~list(X^o, X^a)~']'~" and "~'['~list(X^o, X^a)~']'), "No differences")) +
  theme(text = element_text(size = 18), 
        strip.text = element_text(size = 14),
        legend.box = "vertical",
        legend.key.width = unit(1., "cm"),
        legend.position = "top") + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_boxplot(mapping = aes(ymin = y0, 
                             lower = y25, 
                             middle = y50,
                             upper = y75, 
                             ymax = y100, 
                             fill = method_name,
                             alpha = base_case, 
                             group = interaction(base_case, method_name, n_hist)),
               stat = "identity") + 
  geom_point(mapping = aes(y = y50, group = interaction(base_case, method_name, n_hist), shape = base_case),
             position = position_dodge(width = 0.9)) +
  geom_text(data = rimse_figure_data %>% 
              filter(!different_external_model & different_covariate_dist),
            mapping = aes(label = rank_log10_rel_squared_error_orig,
                          group = method_name, 
                          y = -Inf),
            vjust = -0.2,
            position = position_dodge(width = 0.9)) +
  guides(fill = guide_legend(nrow = 1, order = 0),
         shape = guide_legend(nrow = 1, order = 1),
         alpha = guide_legend(nrow = 1, override.aes = list(fill = "darkgrey"), order = 1))
ggsave(filename = "figures/figure4.png", height = 12, width = 10, dpi = 500)
ggsave(filename = "figures/figure4.eps", height = 12, width = 10)


## -----------------------------

table_mse <- 
  all_scores %>% 
  filter(is_primary | is_secondary) %>%
  filter((!different_external_model&!different_covariate_dist) | ((different_external_model&different_covariate_dist))) %>%
  group_by(scenario_name, true_bivariate_cor, 
           n_curr, n_hist, n_curr_fancy, 
           different_external_model, different_covariate_dist, 
           method_name) %>%
  summarize(mean_log10_rel_squared_error_orig = mean(log10_rel_squared_error_orig, na.rm = TRUE),
            mean_log10_squared_error_orig = mean(log10_squared_error_orig, na.rm = TRUE),
            mean_auc = mean(auc, na.rm = TRUE),
            mean_converge = mean(converge),
            num_datasets = n()) %>%
  mutate(Method =
           case_when(
             method_name == "sab" ~ "SAB", 
             method_name == "sab2" ~ "SAB2",
             method_name == "sab2_alt" ~ "SAB2(Alt)",
             method_name == "sab2_flexible" ~ "SAB2(Ext)",
             method_name == "cml_saddlepoint" ~ "CML(SP)",
             method_name == "cml_newtonraphson" ~ "CML(NR)",
             method_name == "gim_author" ~ "GIM", 
             method_name == "glm_bayes" ~ "GLM(Curr)",
             method_name == "ratios" ~ "Ratios"
           ) %>% factor(levels = c("GLM(Curr)","CML(NR)","CML(SP)", "GIM", "Ratios", "SAB", "SAB2","SAB2(Alt)","SAB2(Ext)")),
         Scen. = 
           case_when(
             scenario_name == "scenario1" ~ "1",
             scenario_name == "scenario2" ~ "2",
             scenario_name == "scenario3" ~ "3",
             scenario_name == "scenario4" ~ "4"
           ) %>% fct_inorder()) %>%
  mutate(logRIMSE = 
           case_when(
             mean_log10_rel_squared_error_orig < min(mean_log10_rel_squared_error_orig) + log(1.05) ~ glue("\\textbf{{{formatC(mean_log10_rel_squared_error_orig,format='f',digits=3)}}}"),
             TRUE ~ formatC(mean_log10_rel_squared_error_orig,format='f',digits=3)),
         logRIMSE = 
           case_when(
             mean_converge >= 1 - .Machine$double.eps ~ logRIMSE,
             mean_converge >= 0.99 ~ paste0(logRIMSE, "*"),
             mean_converge >= 0.96 ~ paste0(logRIMSE, "**"),
             mean_converge >= 0.93 ~ paste0(logRIMSE, "***"),
             TRUE ~ NA_character_
           ),
         logMSE = 
           case_when(
             mean_log10_squared_error_orig < min(mean_log10_squared_error_orig) + log(1.05) ~ glue("\\textbf{{{formatC(mean_log10_squared_error_orig,format='f',digits=3)}}}"),
             TRUE ~ formatC(mean_log10_squared_error_orig,format='f',digits=3)),
         logMSE = 
           case_when(
             mean_converge < 1 ~ paste0(logMSE, "*"), 
             TRUE ~ logMSE
           ),
         #value = glue("{mean_log10_rel_squared_error_orig_pretty}({formatC(100*mean_converge,format='f',digits=0)}\\%)"), 
         #Conv. = glue("{formatC(100*mean_converge,format='f',digits=0)}\\%"),
         AUC = formatC(mean_auc, format = "f", digits = 3),
         different_dist = 
           case_when(
             !different_external_model ~ "same",
             TRUE ~ "different"
           ),
         n_hist = glue("$n_h={n_hist}$") %>% as.character(),
         .keep = "unused") %>%
  ungroup() %>%
  arrange(desc(different_dist), n_curr, desc(n_hist)) 

#### Table 2
#### logRIMSE
table_mse %>%
  pivot_wider(id_cols = c(Scen., Method),
              names_from = c(different_dist,n_curr, n_hist),
              values_from = logRIMSE) %>%
  knitr::kable(format = "latex",
               col.names =  str_remove_all(colnames(.), "[a-z,A-Z]{3,}_[0-9]{3,}_"),
               booktabs = T,
               longtable = F,
               escape = F,
               linesep = c("", "", "", "", "", "", "\\addlinespace")) %>%
  kable_styling(latex_options = c("HOLD_position"),
                full_width = F,
                font_size = 11) %>%
  add_header_above(header = c(" " = 2, "$n=200$" = 2, "$n=400$" = 2, "$n=200$" = 2, "$n=400$" = 2), 
                   escape = F) %>%
  add_header_above(header = c(" " = 2, "Same $[Y|X^o,X^a]$ and $[X^o,X^a]$" = 4, "Different $[Y|X^o,X^a]$ and $[X^o,X^a]$" = 4), 
                   escape = F)

#### Table S2
#### logMSE 

table_mse %>%
  pivot_wider(id_cols = c(Scen., Method),
              names_from = c(different_dist,n_curr, n_hist),
              values_from = logMSE) %>%
  knitr::kable(format = "latex",
               col.names =  str_remove_all(colnames(.), "[a-z,A-Z]{3,}_[0-9]{3,}_"),
               booktabs = T,
               longtable = F,
               escape = F,
               linesep = c("", "", "", "", "", "", "\\addlinespace")) %>%
  kable_styling(latex_options = c("HOLD_position"),
                full_width = F,
                font_size = 11) %>%
  add_header_above(header = c(" " = 2, "$n=200$" = 2, "$n=400$" = 2, "$n=200$" = 2, "$n=400$" = 2), 
                   escape = F) %>%
  add_header_above(header = c(" " = 2, "Same $[Y|X^o,X^a]$ and $[X^o,X^a]$" = 4, "Different $[Y|X^o,X^a]$ and $[X^o,X^a]$" = 4), 
                   escape = F)


## -----------------------------

table_bias_se <- 
  all_betas %>% 
  filter(is_primary | is_secondary) %>%
  filter((!different_external_model&!different_covariate_dist) | ((different_external_model&different_covariate_dist))) %>%
  filter(orig) %>%
  group_by(scenario_name, true_bivariate_cor, 
           n_curr, n_hist, n_curr_fancy, 
           different_external_model, different_covariate_dist, 
           beta_label,
           method_name) %>%
  summarize(rel_bias = 100 * (mean(est_betas, na.rm = T) - first(true_betas)) / abs(first(true_betas)),
            bias = mean(est_betas, na.rm = T) - first(true_betas),
            std_error = sd(est_betas, na.rm = T),
            mean_est_betas = mean(est_betas, na.rm = T),
            mean_converge = mean(converge),
            true_betas = first(true_betas)) %>%
  filter(beta_label == "p1"|
           scenario_name == "scenario2"|
           (scenario_name == "scenario4" & beta_label == "p6")) %>%
  #arrange(desc(abs(true_betas)), beta_label) %>%
  #slice(1) %>% 
  #ungroup(method_name) %>%
  mutate(Method =
           case_when(
             method_name == "sab" ~ "SAB", 
             method_name == "sab2" ~ "SAB2",
             method_name == "sab2_alt" ~ "SAB2(Alt)",
             method_name == "sab2_flexible" ~ "SAB2(Ext)",
             method_name == "cml_saddlepoint" ~ "CML(SP)",
             method_name == "cml_newtonraphson" ~ "CML(NR)",
             method_name == "gim_author" ~ "GIM", 
             method_name == "glm_bayes" ~ "GLM(Curr)",
             method_name == "ratios" ~ "Ratios"
           ) %>% factor(levels = c("GLM(Curr)","CML(NR)","CML(SP)", "GIM", "Ratios", "SAB", "SAB2","SAB2(Alt)","SAB2(Ext)")),
         Scen. = 
           case_when(
             scenario_name == "scenario1" ~ "1",
             scenario_name == "scenario2" ~ "2",
             scenario_name == "scenario3" ~ "3",
             scenario_name == "scenario4" ~ "4"
           ) %>% fct_inorder()) %>%
  mutate(`Rel. Bias` = 
           case_when(
             abs(rel_bias) < min(abs(rel_bias)) + 5 ~ glue("\\textbf{{{formatC(rel_bias,format='f',digits = 1)}}}"),
             TRUE ~ formatC(rel_bias,format='f',digits = 1)),
         `Rel. Bias` = 
           case_when(
             mean_converge < 1 ~ paste0(`Rel. Bias`, "*"), 
             TRUE ~ `Rel. Bias`
           ),
         SE = 
           case_when(
             std_error < 1.05 * min(std_error) ~ glue("\\textbf{{{formatC(std_error,format='f',digits = 3)}}}"),
             TRUE ~ formatC(std_error,format='f',digits = 3)),
         SE = 
           case_when(
             mean_converge < 1 ~ paste0(SE, "*"), 
             TRUE ~ SE
           ),
         different_dist = 
           case_when(
             !different_external_model ~ "same",
             TRUE ~ "different"
           ),
         n_hist = glue("$n_h={n_hist}$") %>% as.character(),
         `$\\beta^o_i$` = true_betas,
         `$i$` = str_remove_all(beta_label, "p"),
         .keep = "unused") %>%
  ungroup() %>%
  arrange(desc(different_dist), n_curr, desc(n_hist)) 

# Table S3
# Relative bias
table_bias_se %>%
  pivot_wider(id_cols = c(Scen., `$i$`, `$\\beta^o_i$`, Method),
              names_from = c(different_dist,n_curr, n_hist),
              values_from = `Rel. Bias`) %>%
  arrange(Scen., desc(`$\\beta^o_i$`)) %>%
  knitr::kable(format = "latex",
               col.names =  str_remove_all(colnames(.), "[a-z,A-Z]{3,}_[0-9]{3,}_"),
               booktabs = T,
               longtable = F,
               escape = F,
               linesep = c("", "", "", "", "", "", "\\addlinespace")) %>%
  kable_styling(latex_options = c("HOLD_position"),
                full_width = F,
                font_size = 11) %>%
  add_header_above(header = c(" " = 4, "$n=200$" = 2, "$n=400$" = 2, "$n=200$" = 2, "$n=400$" = 2), 
                   escape = F) %>%
  add_header_above(header = c(" " = 4, "Same $[Y|X^o,X^a]$ and $[X^o,X^a]$" = 4, "Different $[Y|X^o,X^a]$ and $[X^o,X^a]$" = 4), 
                   escape = F)

# Table S4
# Standard Errors
table_bias_se %>%
  pivot_wider(id_cols = c(Scen., `$i$`, `$\\beta^o_i$`, Method),
              names_from = c(different_dist,n_curr, n_hist),
              values_from = SE) %>%
  arrange(Scen., desc(`$\\beta^o_i$`)) %>%
  knitr::kable(format = "latex",
               col.names =  str_remove_all(colnames(.), "[a-z,A-Z]{3,}_[0-9]{3,}_"),
               booktabs = T,
               longtable = F,
               escape = F,
               linesep = c("", "", "", "", "", "", "\\addlinespace")) %>%
  kable_styling(latex_options = c("HOLD_position"),
                full_width = F,
                font_size = 11) %>%
  add_header_above(header = c(" " = 4, "$n=200$" = 2, "$n=400$" = 2, "$n=200$" = 2, "$n=400$" = 2), 
                   escape = F) %>%
  add_header_above(header = c(" " = 4, "Same $[Y|X^o,X^a]$ and $[X^o,X^a]$" = 4, "Different $[Y|X^o,X^a]$ and $[X^o,X^a]$" = 4), 
                   escape = F)




## -----------------------------
table_conv_running <- 
  all_scores %>% 
  filter(is_primary | is_secondary) %>%
  filter((!different_external_model&!different_covariate_dist) | ((different_external_model&different_covariate_dist))) %>%
  group_by(scenario_name, true_bivariate_cor, 
           n_curr, n_hist, n_curr_fancy, 
           different_external_model, different_covariate_dist, 
           method_name) %>%
  summarize(mean_run_time_secs = mean(run_time),
            sd_run_time_secs = sd(run_time),
            sum_converge = sum(converge),
            num_datasets = n()) %>%
  mutate(Method =
           case_when(
             method_name == "sab" ~ "SAB", 
             method_name == "sab2" ~ "SAB2",
             method_name == "sab2_alt" ~ "SAB2(Alt)",
             method_name == "sab2_flexible" ~ "SAB2(Ext)",
             method_name == "cml_saddlepoint" ~ "CML(SP)",
             method_name == "cml_newtonraphson" ~ "CML(NR)",
             method_name == "gim_author" ~ "GIM", 
             method_name == "glm_bayes" ~ "GLM(Curr)",
             method_name == "ratios" ~ "Ratios"
           ) %>% factor(levels = c("GLM(Curr)","CML(NR)","CML(SP)", "GIM", "Ratios", "SAB", "SAB2","SAB2(Alt)","SAB2(Ext)")),
         Scen. = 
           case_when(
             scenario_name == "scenario1" ~ "1",
             scenario_name == "scenario2" ~ "2",
             scenario_name == "scenario3" ~ "3",
             scenario_name == "scenario4" ~ "4"
           ) %>% fct_inorder()) %>%
  mutate(run_time_val = glue("{formatC(mean_run_time_secs, format = 'f', digits = 0)}({formatC(sd_run_time_secs, format = 'f', digits = 0)})"),
         different_dist = 
           case_when(
             !different_external_model ~ "same",
             TRUE ~ "different"
           ),
         n_hist = glue("$n_h={n_hist}$") %>% as.character(),
         .keep = "unused") %>%
  ungroup() %>%
  arrange(desc(different_dist), n_curr, desc(n_hist)) 

# Table S5
# Convergence
table_conv_running %>%
  pivot_wider(id_cols = c(Scen., Method),
              names_from = c(different_dist,n_curr, n_hist),
              values_from = sum_converge) %>%
  knitr::kable(format = "latex",
               col.names =  str_remove_all(colnames(.), "[a-z,A-Z]{3,}_[0-9]{3,}_"),
               booktabs = T,
               longtable = F,
               escape = F,
               linesep = c("", "", "", "", "", "", "\\addlinespace")) %>%
  kable_styling(latex_options = c("HOLD_position"),
                full_width = F,
                font_size = 11) %>%
  add_header_above(header = c(" " = 2, "$n=200$" = 2, "$n=400$" = 2, "$n=200$" = 2, "$n=400$" = 2), 
                   escape = F) %>%
  add_header_above(header = c(" " = 2, "Same $[Y|X^o,X^a]$ and $[X^o,X^a]$" = 4, "Different $[Y|X^o,X^a]$ and $[X^o,X^a]$" = 4), 
                   escape = F)


# Table S6
# Run time
table_conv_running %>%
  pivot_wider(id_cols = c(Scen., Method),
              names_from = c(different_dist,n_curr, n_hist),
              values_from = run_time_val) %>%
  knitr::kable(format = "latex",
               col.names =  str_remove_all(colnames(.), "[a-z,A-Z]{3,}_[0-9]{3,}_"),
               booktabs = T,
               longtable = F,
               escape = F,
               linesep = c("", "", "", "", "", "", "\\addlinespace")) %>%
  kable_styling(latex_options = c("HOLD_position"),
                full_width = F,
                font_size = 11) %>%
  add_header_above(header = c(" " = 2, "$n=200$" = 2, "$n=400$" = 2, "$n=200$" = 2, "$n=400$" = 2), 
                   escape = F) %>%
  add_header_above(header = c(" " = 2, "Same $[Y|X^o,X^a]$ and $[X^o,X^a]$" = 4, "Different $[Y|X^o,X^a]$ and $[X^o,X^a]$" = 4), 
                   escape = F)



