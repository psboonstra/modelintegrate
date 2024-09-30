
fxn_genmeta_logistic_author <-
  function(X,
           beta_internal_with_intercept,
           theta_tilde_with_intercept, 
           n_hist, 
           tol, max_rep = 1000) {
    
    stopifnot(all(colnames(X)==names(beta_internal_with_intercept)[-1]))
    stopifnot(all(names(theta_tilde_with_intercept) %in% names(beta_internal_with_intercept)))
    
    studies <- list(list(Coeff = beta_internal_with_intercept,
                         Covariance = NULL,
                         Sample_size = nrow(X)), 
                    list(Coeff = theta_tilde_with_intercept,
                         Covariance = NULL,
                         Sample_size = n_hist))
    
    foo <- GENMETA_phil(study_info = studies, 
                        ref_dat = cbind("(Intercept)" = 1, X), 
                        model = "logistic",
                        variable_intercepts = TRUE, 
                        control = list(epsilon = tol, maxit = max_rep)) %>%
      try(silent = TRUE)
    
    if("try-error" %in% class(foo)) {
      
      return(list(intercept_hat = NA,
                  beta_hat = NA,
                  theta_tilde_with_intercept = theta_tilde_with_intercept,
                  converge = FALSE, 
                  message = foo[[1]],
                  iter = NA,
                  singular_hessian = NA,
                  final_diff = NA,
                  objective_function = c("total" = NA))) 
    } else if(any(is.na(foo$Est.coeff))) {
      
      return(list(intercept_hat = NA,
                  beta_hat = NA,
                  theta_tilde_with_intercept = theta_tilde_with_intercept,
                  converge = FALSE, 
                  iter = foo$outer_iter,
                  singular_hessian = NA,
                  final_diff = foo$eps_outer,
                  objective_function = c("total" = NA)))
    } else {
      
      return(list(intercept_hat = foo$Est.coeff["(Intercept_Study_1)"] %>% `names<-`("(Intercept)"),
                  beta_hat = foo$Est.coeff[colnames(X)],
                  theta_tilde_with_intercept = theta_tilde_with_intercept,
                  converge = (foo$outer_iter <= max_rep), 
                  iter = foo$outer_iter,
                  singular_hessian = NA,
                  final_diff = foo$eps_outer,
                  objective_function = c("total" = NA)))
      
    }
  }

# GENMETA_phil is nearly identical to GENMETA::GENMETA, the only changes being
# (i) that there is an added check to see if the number of iterations has exceeded
# the value of maxit and to break if so and (ii) to additionally return the last value
# of eps_outer

GENMETA_phil = function (study_info, ref_dat, model, variable_intercepts = FALSE, 
                         initial_val = NULL, control = list(epsilon = 1e-06, maxit = 1000)) 
{
  call_MetaG <- match.call()
  if (!missing(control)) {
    control = control
  }
  threshold <- control[[1]]
  maxit <- control[[2]]
  different_intercept <- variable_intercepts
  study_estimates <- study_info
  error_1 <- 0
  error_2 <- 0
  error_3 <- 0
  error_4 <- 0
  error_5 <- 0
  no_of_studies <- length(study_estimates)
  temp <- c()
  indicator_missing_covariance_sample_size = 0
  missing_study_sample_size <- c()
  missing_covariance_study_indices <- c()
  for (i in 1:no_of_studies) {
    if (is.null(study_estimates[[i]][[3]]) == T) 
      missing_study_sample_size = c(missing_study_sample_size, 
                                    i)
  }
  for (i in 1:no_of_studies) {
    temp <- union(temp, names(study_estimates[[i]][[1]]))
  }
  for (i in 1:no_of_studies) {
    if (is.null(study_estimates[[i]][[2]]) == T && is.null(study_estimates[[i]][[3]]) == 
        T) 
      indicator_missing_covariance_sample_size = 1
  }
  if (indicator_missing_covariance_sample_size == 1) {
    print("Error: All the studies should have either an estimate for the var-cov matrix or the sample size. Atleast one of them is missing(NULL) in atleast one of the studies")
    error_1 <- 1
  }
  if (indicator_missing_covariance_sample_size == 0) {
    for (i in missing_study_sample_size) {
      if (length(study_estimates[[i]][[1]]) != ncol(study_estimates[[i]][[2]])) {
        print("Error: length of the study parameter(effect sizes) vector does not match with the dimension of its variance covariance matrix")
        error_2 <- 1
      }
      if (sum(names(study_estimates[[i]][[1]]) != colnames(study_estimates[[i]][[2]])) > 
          0) {
        print("Error: names of the variables corresponding to the study specific parameter vector is not same(also not in the same order) as in its variance covariance matrix")
        error_2 <- 1
      }
    }
  }
  if (ncol(ref_dat) < length(temp)) {
    print("Error: number of covariates in the reference data does not match with that of the maximal model")
    error_3 <- 1
  }
  if (sum(is.na(match(temp, colnames(ref_dat)))) > 0) {
    print("Error: names of covariates in the reference data does not match with that of the study specific covariates")
    error_4 <- 1
  }
  if (model == "linear" & variable_intercepts == "TRUE") {
    print("Error: when the model is linear, the current version works only when intercepts are assumed same across studies ")
    error_5 <- 1
  }
  names_wo_intercept <- c()
  for (i in 1:no_of_studies) {
    names_wo_intercept <- union(names_wo_intercept, names(study_estimates[[i]][[1]][-1]))
  }
  weight_sum <- 0
  sum <- 0
  estimates_in_which_studies_indices <- list()
  for (k in 1:length(names_wo_intercept)) {
    temp_estimates_in_which <- c()
    for (j in 1:no_of_studies) {
      if (names_wo_intercept[k] %in% names(study_estimates[[j]][[1]]) == 
          T) 
        temp_estimates_in_which <- c(temp_estimates_in_which, 
                                     j)
    }
    estimates_in_which_studies_indices[[k]] <- temp_estimates_in_which
  }
  if (error_1 == 0 && error_2 == 0 && error_3 == 0 && error_4 == 
      0 && error_5 == 0) {
    if (length(initial_val) == 0) {
      initial_val <- c()
      if (different_intercept == TRUE) {
        initial_val <- c(initial_val, unlist(lapply(lapply(study_estimates, 
                                                           `[[`, 1), `[[`, 1)))
        for (k in 1:length(names_wo_intercept)) {
          for (j in estimates_in_which_studies_indices[[k]]) {
            if (is.null(study_estimates[[j]][[2]]) == 
                F) {
              index_cov <- which(names(study_estimates[[j]][[1]]) %in% 
                                   names_wo_intercept[k] == T)
              weight <- 1/study_estimates[[j]][[2]][index_cov, 
                                                    index_cov]
              weight_sum <- weight_sum + weight
              sum <- sum + study_estimates[[j]][[1]][index_cov] * 
                weight
            }
            if (is.null(study_estimates[[j]][[2]]) == 
                T) {
              index_cov <- which(names(study_estimates[[j]][[1]]) %in% 
                                   names_wo_intercept[k] == T)
              weight <- study_estimates[[j]][[3]]
              weight_sum <- weight_sum + weight
              sum <- sum + study_estimates[[j]][[1]][index_cov] * 
                weight
            }
          }
          initial_val[k + length(study_estimates)] <- sum/weight_sum
          weight_sum <- 0
          sum <- 0
        }
      }
      if (different_intercept == F) {
        for (k in 1:length(names_wo_intercept)) {
          for (j in estimates_in_which_studies_indices[[k]]) {
            if (is.null(study_estimates[[j]][[2]]) == 
                F) {
              index_cov <- which(names(study_estimates[[j]][[1]]) %in% 
                                   names_wo_intercept[k] == T)
              weight <- 1/study_estimates[[j]][[2]][index_cov, 
                                                    index_cov]
              weight_sum <- weight_sum + weight
              sum <- sum + study_estimates[[j]][[1]][index_cov] * 
                weight
            }
            if (is.null(study_estimates[[j]][[2]]) == 
                T) {
              index_cov <- which(names(study_estimates[[j]][[1]]) %in% 
                                   names_wo_intercept[k] == T)
              weight <- study_estimates[[j]][[3]]
              weight_sum <- weight_sum + weight
              sum <- sum + study_estimates[[j]][[1]][index_cov] * 
                weight
            }
          }
          initial_val[k + 1] <- sum/weight_sum
          weight_sum <- 0
          sum <- 0
        }
        for (j in 1:no_of_studies) {
          if (is.null(study_estimates[[j]][[2]]) == F) {
            weight <- 1/study_estimates[[j]][[2]][1, 
                                                  1]
            weight_sum <- weight_sum + weight
            sum <- sum + study_estimates[[j]][[1]][1] * 
              weight
          }
          if (is.null(study_estimates[[j]][[2]]) == T) {
            weight <- study_estimates[[j]][[3]]
            weight_sum <- weight_sum + weight
            sum <- sum + study_estimates[[j]][[1]][1] * 
              weight
          }
        }
        initial_val[1] <- sum/weight_sum
      }
    }
    for (i in 1:no_of_studies) {
      if (is.null(study_estimates[[i]][[2]]) == T) 
        missing_covariance_study_indices = c(missing_covariance_study_indices, 
                                             i)
    }
    X_rbind <- c()
    for (k in 1:no_of_studies) {
      X_rbind <- rbind(X_rbind, ref_dat)
    }
    X_bdiag_list <- list()
    for (k in 1:no_of_studies) {
      col_ind <- which(colnames(ref_dat) %in% names(study_estimates[[k]][[1]]) == 
                         TRUE)
      X_bdiag_list[[k]] <- as.matrix(ref_dat[, col_ind])
    }
    X_bdiag <- Reduce(magic::adiag, X_bdiag_list)
    model_optim <- model
    study_optim <- study_estimates
    eps_outer = 0
    if (model_optim == "linear") {
      study_indices <- seq(1, no_of_studies, 1)
      if (length(missing_covariance_study_indices) > 0) 
        non_missing_covariance_study_indices <- study_indices[-which(study_indices %in% 
                                                                       missing_covariance_study_indices)]
      if (length(missing_covariance_study_indices) == 0) 
        non_missing_covariance_study_indices = study_indices
      disp <- rep(NA, no_of_studies)
      for (j in study_indices) {
        col_ind <- which(colnames(ref_dat) %in% names(study_estimates[[j]][[1]]) == 
                           TRUE)
        disp[j] <- (1 - ((nrow(ref_dat) - 1)/nrow(ref_dat) * 
                           var(ref_dat[, col_ind] %*% study_estimates[[j]][[1]])))
      }
      C_init = diag(ncol(X_bdiag))
      Gamma_hat <- matrix(NA, nrow(C_init), ncol(ref_dat))
      lambda_ref <- list()
      k_gamma_hat <- 1
      for (k in 1:no_of_studies) {
        col_ind <- which(colnames(ref_dat) %in% names(study_estimates[[k]][[1]]) == 
                           TRUE)
        Gamma_hat[k_gamma_hat:(k_gamma_hat + length(col_ind) - 
                                 1), ] <- (t(ref_dat[, col_ind]) %*% ref_dat)/(disp[[k]] * 
                                                                                 nrow(ref_dat))
        k_gamma_hat <- k_gamma_hat + length(col_ind)
      }
      for (k in non_missing_covariance_study_indices) {
        col_ind <- which(colnames(ref_dat) %in% names(study_estimates[[k]][[1]]) == 
                           TRUE)
        temp_lambda_ref <- t(ref_dat[, col_ind]) %*% 
          ref_dat[, col_ind]
        lambda_ref[[k]] <- (temp_lambda_ref %*% study_estimates[[k]][[2]] %*% 
                              temp_lambda_ref)/(nrow(ref_dat) * disp[k]^2)
      }
      A_n1 <- matrix(NA, nrow(C_init), ncol(ref_dat))
      B_n1 <- matrix(NA, nrow(C_init), 1)
      k_A = 1
      for (k in 1:no_of_studies) {
        col_ind <- which(colnames(ref_dat) %in% names(study_estimates[[k]][[1]]) == 
                           TRUE)
        A_n1[k_A:(k_A + length(col_ind) - 1), ] <- (t(ref_dat[, 
                                                              col_ind]) %*% ref_dat) * (1/disp[k])
        B_n1[k_A:(k_A + length(col_ind) - 1), ] <- (t(ref_dat[, 
                                                              col_ind]) %*% ref_dat[, col_ind] %*% study_estimates[[k]][[1]]) * 
          (1/disp[k])
        k_A <- k_A + length(col_ind)
      }
      beta_old_first <- solve(t(A_n1) %*% C_init %*% A_n1, 
                              tol = 1e-60)
      beta_old_sec <- t(A_n1) %*% C_init %*% B_n1
      beta_old_identity <- beta_old_first %*% beta_old_sec
      beta_old <- beta_old_first %*% beta_old_sec
      A_n2 <- 1/(disp^2)
      B_n2 <- rep(NA, no_of_studies)
      proceed <- TRUE
      U <- matrix(NA, nrow(C_init), nrow(ref_dat))
      no_of_iter <- 0
      while (proceed) {
        no_of_iter <- no_of_iter + 1
        disp_max_old <- (1 - ((nrow(ref_dat) - 1)/nrow(ref_dat) * 
                                var(ref_dat %*% beta_old)))
        for (k in missing_covariance_study_indices) {
          col_ind <- which(colnames(ref_dat) %in% names(study_estimates[[k]][[1]]) == 
                             TRUE)
          temp_lambda_ref_1 <- as.numeric(disp_max_old) * 
            t(ref_dat[, col_ind]) %*% ref_dat[, col_ind]
          temp_lambda_W <- diag((as.vector(ref_dat %*% 
                                             beta_old) - as.vector(ref_dat[, col_ind] %*% 
                                                                     study_estimates[[k]][[1]])))
          temp_lambda_ref_2 <- t(ref_dat[, col_ind]) %*% 
            temp_lambda_W %*% ref_dat[, col_ind]
          lambda_ref[[k]] <- (temp_lambda_ref_1 + temp_lambda_ref_2)/(study_estimates[[k]][[3]] * 
                                                                        (disp[k]^2))
        }
        Lambda_ref <- Reduce(magic::adiag, lambda_ref)
        k_U <- 1
        for (k in 1:no_of_studies) {
          col_ind <- which(colnames(ref_dat) %in% names(study_estimates[[k]][[1]]) == 
                             TRUE)
          temp_delta_k <- diag(as.vector(ref_dat %*% 
                                           beta_old) - as.vector(ref_dat[, col_ind] %*% 
                                                                   study_estimates[[k]][[1]]))/(disp[k]^2)
          U[k_U:(k_U + length(col_ind) - 1), ] <- t(ref_dat[, 
                                                            col_ind]) %*% temp_delta_k
          k_U <- k_U + length(col_ind)
        }
        Delta_hat <- (U %*% t(U))/(nrow(ref_dat))
        C_new <- solve(Lambda_ref + Delta_hat, tol = 1e-60)
        asy_var_C_identity <- (solve(t(Gamma_hat) %*% 
                                       Gamma_hat, tol = 1e-60) %*% (t(Gamma_hat) %*% 
                                                                      C_new %*% Gamma_hat) %*% solve(t(Gamma_hat) %*% 
                                                                                                       Gamma_hat, tol = 1e-60))/(nrow(ref_dat))
        beta_new_first <- solve(t(A_n1) %*% C_new %*% 
                                  A_n1, tol = 1e-60)
        beta_new_sec <- t(A_n1) %*% C_new %*% B_n1
        beta_new <- beta_new_first %*% beta_new_sec
        eps_outer = sqrt(sum((beta_new - beta_old)^2))
        beta_old <- beta_new
        if (eps_outer < threshold || no_of_iter > maxit) 
          proceed <- FALSE
      }
      asy_var_opt <- solve(t(Gamma_hat) %*% C_new %*% Gamma_hat, 
                           tol = 1e-60)/nrow(ref_dat)
      beta_old_identity <- as.vector(beta_old_identity)
      names(beta_old_identity) <- colnames(ref_dat)
      if (is.null(asy_var_C_identity) == FALSE) {
        colnames(asy_var_C_identity) <- colnames(ref_dat)
        rownames(asy_var_C_identity) <- colnames(ref_dat)
      }
      if (is.null(asy_var_opt) == FALSE) {
        colnames(asy_var_opt) <- colnames(ref_dat)
        rownames(asy_var_opt) <- colnames(ref_dat)
      }
      beta_old <- as.vector(beta_old)
      names(beta_old) <- colnames(ref_dat)
      linear_result <- list(Est.coeff = beta_old,
                            Est.var.cov = asy_var_opt, 
                            Res.var = disp_max_old, 
                            iter = no_of_iter, 
                            outer_iter = no_of_iter,
                            eps_outer = eps_outer, 
                            call = call_MetaG)
      class(linear_result) <- "GENMETA"
      return(linear_result)
    }
    if (model_optim == "logistic") {
      C_init = diag(ncol(X_bdiag))
      no_of_iter_outer = 0
      total_iter = 0
      output_identity <- GENMETA:::myoptim(no_of_studies, study_optim, 
                                           ref_dat, X_rbind, X_bdiag_list, C_init, initial_val, 
                                           threshold, model_optim, missing_covariance_study_indices, 
                                           different_intercept, no_of_iter_outer)
      beta_identity <- output_identity$beta_optim
      beta_initial <- output_identity$beta_optim
      asy_var_beta_converged_identity <- output_identity$Asy_var_optim
      total_iter_identity <- output_identity$iter_IRWLS
      C_iter <- output_identity$C_optim
      proceed <- TRUE
      mark_failure = 0
      if (sum(is.na(beta_identity)) > 0 || output_identity$Status == 
          0) {
        beta_initial = rep(NA, ncol(ref_dat))
        asy_var_beta_converged = NULL
      }
      else {
        while (proceed) {
          no_of_iter_outer <- no_of_iter_outer + 1
          output_optim <- GENMETA:::myoptim(no_of_studies, study_optim, 
                                            ref_dat, X_rbind, X_bdiag_list, C_iter, beta_initial, 
                                            threshold, model_optim, missing_covariance_study_indices, 
                                            different_intercept, no_of_iter_outer)
          beta_iter_old <- output_optim$beta_optim
          C_iter <- output_optim$C_optim
          if (sum(is.na(beta_iter_old)) > 0 || output_optim$Status == 
              0) {
            mark_failure = 1
            break
          }
          eps_outer = sqrt(sum((beta_iter_old - beta_initial)^2))
          total_iter <- output_optim$iter_IRWLS + no_of_iter_outer
          beta_initial <- beta_iter_old
          if (eps_outer < threshold || no_of_iter_outer > maxit) 
            proceed <- FALSE
        }
        total_iter <- total_iter + total_iter_identity
        if (mark_failure == 1) {
          beta_initial = rep(NA, ncol(ref_dat))
          asy_var_beta_converged <- NULL
        }
        else {
          asy_var_beta_converged <- output_optim$Asy_var_optim
          if (different_intercept == T) {
            temp_name <- c()
            for (i in 1:no_of_studies) {
              temp_name <- c(temp_name, paste0("(Intercept_Study_", 
                                               i, ")"))
            }
            temp_name <- c(temp_name, colnames(ref_dat[, 
                                                       -1]))
            beta_identity <- as.vector(beta_identity)
            names(beta_identity) <- temp_name
            beta_initial <- as.vector(beta_initial)
            names(beta_initial) <- temp_name
            if (is.null(asy_var_beta_converged_identity) == 
                FALSE) {
              colnames(asy_var_beta_converged_identity) <- temp_name
              rownames(asy_var_beta_converged_identity) <- temp_name
            }
            if (is.null(asy_var_beta_converged) == FALSE) {
              colnames(asy_var_beta_converged) <- temp_name
              rownames(asy_var_beta_converged) <- temp_name
            }
          }
          if (different_intercept == F) {
            beta_identity <- as.vector(beta_identity)
            names(beta_identity) <- colnames(ref_dat)
            beta_initial <- as.vector(beta_initial)
            names(beta_initial) <- colnames(ref_dat)
            if (is.null(asy_var_beta_converged_identity) == 
                FALSE) {
              colnames(asy_var_beta_converged_identity) <- colnames(ref_dat)
              rownames(asy_var_beta_converged_identity) <- colnames(ref_dat)
            }
            if (is.null(asy_var_beta_converged) == FALSE) {
              colnames(asy_var_beta_converged) <- colnames(ref_dat)
              rownames(asy_var_beta_converged) <- colnames(ref_dat)
            }
          }
        }
      }
      logistic_result <- list(Est.coeff = beta_initial, 
                              Est.var.cov = asy_var_beta_converged, 
                              Res.var = NA, 
                              iter = total_iter,
                              outer_iter = no_of_iter_outer,
                              eps_outer = eps_outer, 
                              call = call_MetaG)
      class(logistic_result) <- "GENMETA"
      return(logistic_result)
    }
  }
}