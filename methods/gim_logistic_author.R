
fxn_gim_logistic_author <-
  function(Y, 
           X,
           theta_tilde_with_intercept, 
           n_hist) {
    
    stopifnot(all(names(theta_tilde_with_intercept)[-1] %in% colnames(X)))
    
    orig_var_names <- setdiff(names(theta_tilde_with_intercept), "(Intercept)")
    
    model <- list(list(form = glue('Y~{glue_collapse(orig_var_names, sep = "+")}'), 
                       info = data.frame(var = names(theta_tilde_with_intercept), 
                                         bet = as.numeric(theta_tilde_with_intercept))))
    
    foo <- 
      bind_cols(Y = Y, X) %>%
      gim::gim(formula = glue('Y~{glue_collapse(colnames(X), sep = "+")}'),
               family = "binomial", model = model, nsample = matrix(n_hist)) %>%
      try(silent = TRUE)
    
    if("try-error" %in% class(foo)) {
      
      return(list(intercept_hat = NA,
                  beta_hat = NA,
                  theta_tilde_with_intercept = theta_tilde_with_intercept,
                  Sigma_h_inv = NA,
                  converge = FALSE, 
                  message = foo[[1]],
                  iter = NA,
                  singular_hessian = NA,
                  final_diff = NA,
                  objective_function = c("total" = NA))) 
    } else {
      return(list(intercept_hat = foo$coefficients["(Intercept)"],
                  beta_hat = foo$coefficients[colnames(X)],
                  theta_tilde_with_intercept = theta_tilde_with_intercept,
                  Sigma_h_inv = solve(foo$V.bet) / n_hist,
                  converge = TRUE, 
                  message = NA,
                  iter = NA,
                  singular_hessian = NA,
                  final_diff = NA,
                  objective_function = c("total" = NA))) 
    }
    
  }