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
library(logistf);
library(knitr);
library(kableExtra);


## --------------------------------------------------------------
sims_per_scenario = 200;
which_batch = 1;
total_num_arrays_by_batch = c(999, 999, 999)
total_num_arrays = sum(total_num_arrays_by_batch)


## --------------------------------------------------------------
source("sim_functions/OR_from_AUC.R")
source("sim_functions/generate_params.R")
source("sim_functions/scenario_setup.R")
source("sim_functions/score_method.R")
source("aux_functions/log1plex.R")
source("aux_functions/get_num_eff.R")


## --------------------------------------------------------------
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
  iter = col_integer(), 
  converge = col_logical(), 
  singular = col_logical(), 
  final_diff = col_double(),
  final_objective_function = col_double()
  
);

col_types_coefs <- cols(
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
  beta_label = col_integer(),
  est_betas = col_double()
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
  scenario_complexity_weight = col_double(),
  num_arrays = col_double(),
  start_array_id = col_double(),
  end_array_id = col_double(),
  array_id = col_double(),
  n_sim_this_array_id = col_double(),
  total_runtime_secs = col_double()
);


## --------------------------------------------------------------
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
    #raw_all_scores = 
    #  bind_rows(raw_all_scores, foo);
    assign(paste0("job",i,"_scores"), foo)
  } else {
    cat("sim ", i, ", not found\n");
    not_found = c(not_found,i)
  }
  foo <- try(read_csv(paste0("out/job",i,"_coefs.csv"), col_types = col_types_coefs));
  if(!"try-error" %in% class(foo)) {
    #raw_all_coefs = 
    #  bind_rows(raw_all_coefs, foo);
    assign(paste0("job",i,"_coefs"), foo)
  } 
  foo <- try(read_csv(paste0("out/job",i,"_array_stats.csv"), col_types = col_types_array_stats));
  if(!"try-error" %in% class(foo)) {
    #raw_all_coefs = 
    #  bind_rows(raw_all_coefs, foo);
    assign(paste0("job",i,"_array_stats"), foo)
  } 
  
  if(i%%500 == 0) {cat(i, "\n")}
}

rm(col_types_coefs, col_types_scores, col_types_array_stats)

raw_all_scores <-
  bind_rows(map(str_subset(ls(), pattern = "_scores"), get)) %>%
  arrange(array_id, sim_num)
raw_all_coefs <-
  bind_rows(map(str_subset(ls(), pattern = "_coefs"), get)) %>%
  arrange(array_id, sim_num)
raw_all_array_stats <-
  bind_rows(map(str_subset(ls(), pattern = "_array_stats"), get)) %>%
  arrange(array_id, sim_num)
unique(raw_all_scores$method_name)


## --------------------------------------------------------------
raw_all_array_stats %>% 
  group_by(n_hist, n_curr, different_external_model, different_covariate_dist,scenario_name, scenario, start_array_id, end_array_id) %>% 
  summarize(n_sim = sum(n_sim_this_array_id), total_runtime_secs = sum(total_runtime_secs)) %>% 
  left_join(all_scenarios %>% select(scenario, num_orig, num_aug, scenario_complexity_weight, old_scenario_complexity_weight)) %>%
  ggplot() + 
  geom_point(aes(x = total_runtime_secs, 
                 y = scenario_complexity_weight, 
                 color = scenario_name, 
                 shape = factor(different_covariate_dist)))

raw_all_array_stats %>% 
  group_by(n_hist, n_curr, different_external_model, different_covariate_dist,scenario_name, scenario, start_array_id, end_array_id) %>% 
  summarize(n_sim = sum(n_sim_this_array_id), total_runtime_secs = sum(total_runtime_secs)) %>% 
  left_join(all_scenarios %>% select(scenario, num_orig, num_aug, scenario_complexity_weight)) %>%
  lm(formula = log(total_runtime_secs) ~ factor(scenario_name) + log(n_hist) * log(n_curr) ) %>%
  predict()


## --------------------------------------------------------------
primary_method_names = c(
  "cml_saddlepoint", 
  "gim_orig", 
  "ratios", 
  "sab", # eta scales variance
  "sab2"#, # omega scales mean, omega^2 scales variance, linear projection matrix
)

secondary_method_names = c(
  "cml_newtonraphson", 
  "sab2_flexible" # like sab2, but uses internal covariance matrix 
)


## --------------------------------------------------------------
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


## --------------------------------------------------------------
all_scores <- 
  full_join(
    raw_all_scores %>%
      mutate(converge =
               case_when(
                 str_detect(method_name, "bayes") ~ TRUE, 
                 str_detect(method_name, "sab") ~ TRUE, 
                 str_detect(method_name, "truth") ~ TRUE, 
                 TRUE ~ converge
               )),
    raw_all_scores %>%
      filter(method_name == "glm_bayes") %>% 
      rename_with(~ paste0(.x, "_glm_bayes"), c(squared_error, squared_error_orig, squared_error_aug, converge)) %>%
      dplyr::select(array_id:scenario_name, squared_error_glm_bayes: squared_error_aug_glm_bayes, converge_glm_bayes)) %>%
  mutate(rel_squared_error = squared_error / squared_error_glm_bayes, 
         log10_rel_squared_error = log10(rel_squared_error),
         is_primary = method_name %in% primary_method_names, 
         is_secondary = method_name %in% secondary_method_names) %>%
  # whether to group by n_hist in this next line is important, as it dictates
  # whether results that only vary by n_hist should be all included or all excluded
  group_by(n_hist, n_curr, true_bivariate_cor, scenario_name, 
           different_external_model, different_covariate_dist, data_seed) %>%
  #group_by(n_curr, true_bivariate_cor, scenario_name, 
  #         different_external_model, different_covariate_dist, data_seed) %>%
  mutate(all_primary_converge = all(converge | !is_primary)) %>%
  mutate(all_primary_secondary_converge = all(converge | (!is_primary & !is_secondary))) %>%
  ungroup() %>%
  mutate(n_curr_fancy = glue("n=={n_curr}")) %>%
  mutate(n_hist = factor(n_hist))

all_coefs <- 
  raw_all_coefs %>%
  mutate(converge =
           case_when(
             str_detect(method_name, "bayes") ~ TRUE, 
             str_detect(method_name, "sab") ~ TRUE, 
             str_detect(method_name, "truth") ~ TRUE, 
             TRUE ~ converge
           ),
         is_primary = method_name %in% primary_method_names,
         is_secondary = method_name %in% secondary_method_names) %>%
  # whether to group by n_hist in this next line is important:
  # - if you **do** group by n_hist, then each value of n_hist will potentially 
  # have a different set of data_seeds for which all methods converge. Do this
  # if you want to drop fewer data_seeds
  # - if you **don't** group by n_hist, then you will end up with a common set of
  # data_seeds across all vlaues of n_hist. Do this if you have one or more methods
  # the peformance of which does not depend upon n_hist. 
  group_by(n_hist, n_curr, true_bivariate_cor, scenario_name, 
           different_external_model, different_covariate_dist, data_seed) %>%
  # group_by(n_curr, true_bivariate_cor, scenario_name, 
  #         different_external_model, different_covariate_dist, data_seed) %>%
  mutate(all_primary_converge = all(converge | !is_primary)) %>%
  mutate(all_primary_secondary_converge = all(converge | (!is_primary & !is_secondary))) %>%
  ungroup() %>%
  mutate(n_curr_fancy = glue("n=={n_curr}")) %>%
  mutate(n_hist = factor(n_hist))


rimse_figure_data <- 
  all_scores %>% 
  filter(all_primary_converge, 
         is_primary) %>%
  group_by(scenario_name, true_bivariate_cor, 
           n_curr, n_hist, n_curr_fancy, 
           method_name,
           different_external_model, different_covariate_dist) %>%
  mutate(median_log10_rel_squared_error = median(log10_rel_squared_error), 
         num_datasets = n()) %>%
  ungroup(method_name) %>%
  mutate(rank_log10_rel_squared_error = rank(median_log10_rel_squared_error, ties = "max") / num_datasets) %>%
  ungroup() %>%
  mutate(num_datasets = glue("{num_datasets} / {max(num_datasets)}") %>% as.character(),
         method_name =
           case_when(
             method_name == "sab" ~ "SAB", 
             method_name == "sab2" ~ "SAB2",
             method_name == "sab2_alt" ~ "SAB2(Alt)",
             method_name == "sab2_flexible" ~ "SAB2(Ext)",
             method_name == "cml_saddlepoint" ~ "CML(SP)",
             method_name == "cml_newtonraphson" ~ "CML(NR)",
             method_name == "gim_orig" ~ "GIM", 
             method_name == "glm_bayes" ~ "GLM(Curr)",
             method_name == "ratios" ~ "Ratios"
           ) %>% factor(levels = c("GLM(Curr)","CML(NR)","CML(SP)", "GIM", "Ratios", "SAB", "SAB2","SAB2(Alt)","SAB2(Ext)")),
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
             scenario_name == "scenario1" ~ "list(beta^o==group('{',list(0.637,0.637),'}'),beta^a==1.274)",
             scenario_name == "scenario2" ~ "list(beta^o==group('{',list(0.771, 0.309, 0.031, 0.031),'}'),beta^a==group('{',list(-0.463, 0.463),'}'))",
             scenario_name == "scenario3" ~ "list(group('{',list(p,q),'}')==group('{',list(20,10),'}'), uniform~signal)",
             scenario_name == "scenario4" ~ "list(group('{',list(p,q),'}')==group('{',list(30,10),'}'), nonuniform~signal)"
           ) %>% fct_inorder(),
         scenario_name_fancy_simple = 
           case_when(
             scenario_name == "scenario1" ~ glue("Scen~1~':'~group('{{',list(p,q),'}}')==group('{{',list(2,1),'}}')~';'~{n_curr_fancy}"),
             scenario_name == "scenario2" ~ glue("Scen~2~':'~group('{{',list(p,q),'}}')==group('{{',list(4,2),'}}')~';'~{n_curr_fancy}"),
             scenario_name == "scenario3" ~ glue("Scen~3~':'~group('{{',list(p,q),'}}')==group('{{',list(20,10),'}}')~';'~{n_curr_fancy}"),
             scenario_name == "scenario4" ~ glue("Scen~4~':'~group('{{',list(p,q),'}}')==group('{{',list(30,10),'}}')~';'~{n_curr_fancy}")
           ) %>% fct_inorder()) %>%
  group_by(n_hist, n_curr, n_curr_fancy, true_bivariate_cor, different_external_model, different_covariate_dist, 
           scenario_name, scenario_name_fancy_simple, base_case, 
           method_name, 
           num_datasets, rank_log10_rel_squared_error) %>%
  summarize(y0 = quantile(log10_rel_squared_error, 0.25),
            y25 = quantile(log10_rel_squared_error, 0.25), 
            y50 = quantile(log10_rel_squared_error, 0.5), 
            y75 = quantile(log10_rel_squared_error, 0.75), 
            y100 = quantile(log10_rel_squared_error, 0.75)) %>%
  ungroup() 


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
  geom_text(mapping = aes(label = rank_log10_rel_squared_error,
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
                     values = c(1, 0.5),
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
  geom_text(data = rimse_figure_data %>% 
              filter(!different_external_model & different_covariate_dist),
            mapping = aes(label = rank_log10_rel_squared_error,
                          group = method_name, 
                          y = -Inf),
            vjust = -0.2,
            position = position_dodge(width = 0.9)) +
  guides(fill = guide_legend(nrow = 1, order = 0),
         alpha = guide_legend(nrow = 1, override.aes = list(fill = "#1B9E77"), order = 1))
ggsave(filename = "figures/figure2.png", height = 12, width = 10, dpi = 500)



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
                     values = c(1, 0.5),
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
  geom_text(data = rimse_figure_data %>% 
              filter(different_external_model & !different_covariate_dist),
            mapping = aes(label = rank_log10_rel_squared_error,
                          group = method_name, 
                          y = -Inf),
            vjust = -0.2,
            position = position_dodge(width = 0.9)) +
  guides(fill = guide_legend(nrow = 1, order = 0),
         alpha = guide_legend(nrow = 1, override.aes = list(fill = "#1B9E77"), order = 1))
ggsave(filename = "figures/figure3.png", height = 12, width = 10, dpi = 500)



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
                     values = c(1, 0.5),
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
  geom_text(data = rimse_figure_data %>% 
              filter(different_external_model & different_covariate_dist),
            mapping = aes(label = rank_log10_rel_squared_error,
                          group = method_name, 
                          y = -Inf),
            vjust = -0.2,
            position = position_dodge(width = 0.9)) +
  guides(fill = guide_legend(nrow = 1, order = 0),
         alpha = guide_legend(nrow = 1, override.aes = list(fill = "#1B9E77"), order = 1))
ggsave(filename = "figures/figure4.png", height = 12, width = 10, dpi = 500)


table2 <- 
  all_scores %>% 
  filter(is_primary | is_secondary, 
         !different_external_model, !different_covariate_dist) %>%
  group_by(scenario_name, true_bivariate_cor, 
           n_curr, n_hist, n_curr_fancy, 
           method_name) %>%
  summarize(mean_log10_rel_squared_error = mean(log10_rel_squared_error, na.rm = TRUE),
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
             method_name == "gim_orig" ~ "GIM", 
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
  mutate(mean_log10_rel_squared_error_pretty = 
           case_when(
             mean_log10_rel_squared_error < min(mean_log10_rel_squared_error) + log(1.05) ~ glue("\\textbf{{{formatC(mean_log10_rel_squared_error,format='f',digits=3)}}}"),
             TRUE ~ formatC(mean_log10_rel_squared_error,format='f',digits=3)),
         value = glue("{mean_log10_rel_squared_error_pretty}({formatC(100*mean_converge,format='f',digits=0)}\\%)"), 
         n_hist = glue("$n_h={n_hist}$") %>% as.character(),
         .keep = "unused") %>%
  pivot_wider(id_cols = c(Scen., Method),
              names_from = c(n_curr, n_hist)) 

table2 %>%
  knitr::kable(format = "latex",
               col.names =  str_remove(colnames(.), "[0-9]{3,}_"),
               booktabs = T,
               longtable = F,
               escape = F,
               linesep = c("", "", "", "", "", "", "\\addlinespace")) %>%
  kable_styling(latex_options = c("HOLD_position"),
                full_width = F,
                font_size = 11) %>%
  add_header_above(header = c(" " = 2, "$n=100$" = 3, "$n=400$" = 3), 
                   escape = F)


#### Table s1

tableS1 <- 
  all_scores %>% 
  filter(is_primary | is_secondary, 
         different_external_model, different_covariate_dist) %>%
  group_by(scenario_name, true_bivariate_cor, 
           n_curr, n_hist, n_curr_fancy, 
           method_name) %>%
  summarize(mean_log10_rel_squared_error = mean(log10_rel_squared_error, na.rm = TRUE),
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
             method_name == "gim_orig" ~ "GIM", 
             method_name == "glm_bayes" ~ "GLM(Curr)",
             method_name == "ratios" ~ "Ratios"
           ) %>% factor(levels = c("CML(NR)","CML(SP)", "GIM", "Ratios", "SAB", "SAB2","SAB2(Alt)","SAB2(Ext)")),
         Scen. = 
           case_when(
             scenario_name == "scenario1" ~ "1",
             scenario_name == "scenario2" ~ "2",
             scenario_name == "scenario3" ~ "3",
             scenario_name == "scenario4" ~ "4"
           ) %>% fct_inorder()) %>%
  mutate(mean_log10_rel_squared_error_pretty = 
           case_when(
             mean_log10_rel_squared_error < min(mean_log10_rel_squared_error) + log(1.05) ~ glue("\\textbf{{{formatC(mean_log10_rel_squared_error,format='f',digits=3)}}}"),
             TRUE ~ formatC(mean_log10_rel_squared_error,format='f',digits=3)),
         value = glue("{mean_log10_rel_squared_error_pretty}({formatC(100*mean_converge,format='f',digits=0)}\\%)"), 
         n_hist = glue("$n_h={n_hist}$") %>% as.character(),
         .keep = "unused") %>%
  pivot_wider(id_cols = c(Scen., Method),
              names_from = c(n_curr, n_hist)) 

tableS1 %>%
  knitr::kable(format = "latex",
               col.names =  str_remove(colnames(.), "[0-9]{3,}_"),
               booktabs = T,
               longtable = F,
               escape = F,
               linesep = c("", "", "", "", "", "", "\\addlinespace")) %>%
  kable_styling(latex_options = c("HOLD_position"),
                full_width = F,
                font_size = 11) %>%
  add_header_above(header = c(" " = 2, "$n=100$" = 3, "$n=400$" = 3), 
                   escape = F)




