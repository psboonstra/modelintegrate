score_method <- function(fit, true_betas, num_orig, num_aug, x_validate, y_validate) {
  
  if(is.na(fit$intercept_hat) || any(is.na(fit$beta_hat))) {
    return(tibble(squared_error = NA, 
                  squared_error_orig = NA,
                  squared_error_aug = NA, 
                  brier_score = NA,
                  auc = NA,
                  converge = fit$converge,
                  message = fit$message,
                  iter = fit$iter, 
                  singular = fit$singular_hessian,
                  final_diff = fit$final_diff, 
                  final_objective_function = fit$objective_function["total"]))
  } else {
    
    prob_validate = drop(expit(fit$intercept_hat + x_validate%*%fit$beta_hat))
    auc <- as.numeric(pROC::auc(pROC::roc(y_validate ~ prob_validate, levels = c("0","1"), direction = "<")))
    
    return(tibble(squared_error = mean((fit$beta_hat - true_betas)^2), 
                  squared_error_orig = mean((fit$beta_hat - true_betas)[1:num_orig]^2),
                  squared_error_aug = mean((fit$beta_hat - true_betas)[num_orig+ (1:num_aug)]^2), 
                  brier_score = mean((y_validate - prob_validate)^2),
                  auc = auc,
                  converge = fit$converge,
                  message = fit$message,
                  iter = fit$iter, 
                  singular = fit$singular_hessian,
                  final_diff = fit$final_diff, 
                  final_objective_function = fit$objective_function["total"]))
  }
}