

if(!exists("n_hist_seq")) {n_hist_seq = c(400, 1600);}
if(!exists("n_curr_seq")) {n_curr_seq = c(200, 400);}
if(!exists("different_external_model_seq")) {different_external_model_seq = c(F, T);}
if(!exists("different_covariate_dist_seq")) {different_covariate_dist_seq = c(F, T);}

all_scenarios <- 
  expand_grid(n_hist = n_hist_seq, 
              n_curr = n_curr_seq, 
              different_external_model = different_external_model_seq,
              different_covariate_dist = different_covariate_dist_seq,
              scenario_name = names(scenario_list)) %>%
  mutate(
    scenario = 1:n(),
    num_orig = 
      map_int(scenario_list, ~length(.x$orig))[scenario_name], 
    num_aug = 
      map_int(scenario_list, ~length(.x$aug))[scenario_name], 
    old_scenario_complexity_weight = 
      (I(scenario_name == "scenario1") +
         exp(0.18213) * I(scenario_name == "scenario2") +
         exp(2.25265) * I(scenario_name == "scenario3") +
         exp(2.71712) * I(scenario_name == "scenario4")) *
      (n_hist)^(0.51388) * (n_curr)^(0.82694 - 0.05533 * log(n_hist)),
    old_scenario_complexity_weight = old_scenario_complexity_weight / max(old_scenario_complexity_weight),
    scenario_complexity_weight = 
      case_when(
        scenario_name == "scenario1" & n_hist == 400 & n_curr == 200 ~ 40,
        scenario_name == "scenario2" & n_hist == 400 & n_curr == 200 ~ 58,
        scenario_name == "scenario3" & n_hist == 400 & n_curr == 200 ~ 462,
        scenario_name == "scenario4" & n_hist == 400 & n_curr == 200 ~ 802,
        scenario_name == "scenario1" & n_hist == 400 & n_curr == 400 ~ 68,
        scenario_name == "scenario2" & n_hist == 400 & n_curr == 400 ~ 79,
        scenario_name == "scenario3" & n_hist == 400 & n_curr == 400 ~ 663,
        scenario_name == "scenario4" & n_hist == 400 & n_curr == 400 ~ 934,
        scenario_name == "scenario1" & n_hist == 1600 & n_curr == 200 ~ 65,
        scenario_name == "scenario2" & n_hist == 1600 & n_curr == 200 ~ 78,
        scenario_name == "scenario3" & n_hist == 1600 & n_curr == 200 ~ 567,
        scenario_name == "scenario4" & n_hist == 1600 & n_curr == 200 ~ 1039,
        scenario_name == "scenario1" & n_hist == 1600 & n_curr == 400 ~ 100,
        scenario_name == "scenario2" & n_hist == 1600 & n_curr == 400 ~ 100,
        scenario_name == "scenario3" & n_hist == 1600 & n_curr == 400 ~ 809,
        scenario_name == "scenario4" & n_hist == 1600 & n_curr == 400 ~ 1164
      ),
    scenario_complexity_weight = scenario_complexity_weight / max(scenario_complexity_weight),
    # every scenario gets one array to start. the rest are allocated according to perceived complexity
    # more complex scenarios get more arrays
    num_arrays = 1)

if(nrow(all_scenarios) > total_num_arrays_by_batch[which_batch]) {
  stop(glue("The total number of unique scenarios ({nrow(all_scenarios)})
  exceeds the number of arrays to run ({total_num_arrays_by_batch[which_batch]}) but should not. 
  Increase 'total_num_arrays_by_batch[which_batch]' or decrease the number of unique scenarios."))
}
if((sims_per_scenario * nrow(all_scenarios)) < total_num_arrays_by_batch[which_batch]) {
  stop(glue("The total number of unique scenarios ({nrow(all_scenarios)}) 
  times the total number of simulations per scenario ({sims_per_scenario})
  is less than the total number of arrays to run ({total_num_arrays_by_batch[which_batch]}), 
  which means you have unallocated arrays. Increase 'sims_per_scenario' 
  or decrease 'total_num_arrays_by_batch[which_batch]'"))
}

new_num_arrays <- pull(all_scenarios, num_arrays)
unallocated_arrays <- 
  total_num_arrays_by_batch[which_batch] - sum(new_num_arrays)
temp_complexity_weight <-
  pull(all_scenarios, scenario_complexity_weight) * 
  (new_num_arrays < sims_per_scenario)

# setting the seed is important because this script will be called multiple 
# times and we need the same results across instances
set.seed(1)
while(unallocated_arrays > 0 && any(new_num_arrays < sims_per_scenario)) {
  index_to_increment <- 
    table(sample(nrow(all_scenarios), unallocated_arrays, T, temp_complexity_weight))
  new_num_arrays[as.numeric(names(index_to_increment))] = 
    pmin(sims_per_scenario, new_num_arrays[as.numeric(names(index_to_increment))] + index_to_increment)
  unallocated_arrays <- 
    total_num_arrays_by_batch[which_batch] - sum(new_num_arrays)
  temp_complexity_weight <-
    pull(all_scenarios, scenario_complexity_weight) * 
    (new_num_arrays < sims_per_scenario)
}

all_scenarios <- 
  all_scenarios %>%
  mutate(num_arrays = new_num_arrays,
         start_array_id = lag(1 + cumsum(num_arrays), 1, default = 1), 
         end_array_id = cumsum(num_arrays))
rm(unallocated_arrays, new_num_arrays)
rm(n_hist_seq, n_curr_seq, different_external_model_seq, different_covariate_dist_seq)


