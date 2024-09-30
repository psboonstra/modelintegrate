fxn_glm_bayes <- function(Y, X, seed = sample.int(.Machine$integer.max, 1)) {
  
  fitted_model <- stan_glm(Y ~ .,
                           data = data.frame(Y, X),
                           family = "binomial", 
                           prior = student_t(df = 3),
                           warmup = 1e3, 
                           iter = 2e3,
                           chains = 2,
                           cores = parallel::detectCores(),
                           verbose = FALSE, 
                           seed = seed, 
                           refresh = 0);
  
  beta_hat <- coef(fitted_model)[-1]
  sampler_params <- get_sampler_params(fitted_model$stanfit, inc_warmup = F)
  
  
  return(list(
    intercept_hat = coef(fitted_model)["(Intercept)"],
    beta_hat = beta_hat,
    converge = NA, 
    message = NA,
    iter = NA,
    singular_hessian = NA,
    final_diff = NA,
    objective_function = c("total" = NA),
    num_divergences = sum(sapply(sampler_params, function(x) mean(x[, "divergent__"]))),
    max_rhat = max(summary(fitted_model$stanfit)$summary[,"Rhat"])
    ))
  
}