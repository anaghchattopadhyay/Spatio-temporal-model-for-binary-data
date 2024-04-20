for (i in (1:n_testing)) {
    pred_vec = vector()
    for (j in (1:n)) {
      pred_vec = c(pred_vec, exp(-phi_sp *  (euc.dist(
        c(location_mat_training[j, 2],  location_mat_training[j, 1]),
        c(location_mat_testing[i, 2], location_mat_testing[i, 1])
      ))) * exp(-phi_tmp * (time_testing[i] - time_training[j])))#/20
    }
    
    #prediction regression matrix
    prediction_regress_mat = reg_mat[i,]
    
    S12 = matrix(pred_vec, nrow = num_time)
    
    R = matrix(rowMeans(omega_test), nrow = num_time)
    
    pi_pred_test = t(prediction_regress_mat) %*% rowMeans(theta_test) +
      sum(pred_vec * c(inv_tmp %*% R %*% inv_sp)) +  sqrt(mean(sigma_w_test)) *
      t(chol(1 - sum(pred_vec * c(
        inv_tmp %*% S12 %*% inv_sp
      )))) * rnorm(1000) +
      sqrt(mean(Sigma_e_test)) * rnorm(1000)
    
    pi_predicted = mean(as.numeric(pi_pred_test > 0))
    print(pi_predicted)
    print(Categorical_testing[i])
    if (round(pi_predicted) == 1 &&
        Categorical_testing[i] == 1) {
      t1o1 = t1o1 + 1
    }
    if (round(pi_predicted) == 0 &&
        Categorical_testing[i] == 1) {
      t1o0 = t1o0 + 1
    }
    if (round(pi_predicted) == 1 &&
        Categorical_testing[i] == 0) {
      t0o1 = t0o1 + 1
    }
    if (round(pi_predicted) == 0 &&
        Categorical_testing[i] == 0) {
      t0o0 = t0o0 + 1
    }
    print("_")
    pi_abs = abs(pi_predicted - Categorical_testing[i])
    error_abs_vector = c(error_abs_vector, pi_abs)
    pi_p = c(pi_p, pi_predicted)
    
  }
