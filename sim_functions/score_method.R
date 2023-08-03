score_method <- function(fit, true_betas, num_orig, num_aug, x_validate, y_validate) {
  
  tibble(
    squared_error = mean((fit$beta_hat - true_betas)^2), 
    squared_error_orig = mean((fit$beta_hat - true_betas)[1:num_orig]^2),
    squared_error_aug = mean((fit$beta_hat - true_betas)[num_orig+ (1:num_aug)]^2), 
    brier_score = mean((y_validate - expit(fit$intercept_hat + x_validate%*%fit$beta_hat))^2),
    iter = fit$iter, 
    converge = fit$converge,
    singular = fit$singular_hessian,
    final_diff = fit$final_diff, 
    final_objective_function = fit$objective_function["total"])
  
}