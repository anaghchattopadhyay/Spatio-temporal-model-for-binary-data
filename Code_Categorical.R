  data_modified = training
  location_mat = cbind(data_modified$lat, data_modified$long)
  #have columns as "lat" and "long" corr to training data
  colnames(location_mat) = c("lat", "long")
  
  key = data_modified$id
  
  phi_sp = 1
  phi_tmp = 0.1
  data_modified_size = 10000
  
  gap = 10
  
  GR_threshold = 1.2
  
  max_burn = 100000
  
  inv_gamma_parameter_1 = 2
  
  inv_gamma_parameter_2 = 1
  
  time_points_common_used = c(1:((
    as.numeric(as.POSIXct(test_start[iter_num])) - min(as.numeric(as.POSIXct(dff$week)))
  ) / (86400 * 7)))
  #timepoints common used for training
  tttime_points_common_used = time_points_common_used = c(1:((
    as.numeric(as.POSIXct(test_start[iter_num])) - min(as.numeric(as.POSIXct(dff$week)))
  ) / (86400 * 7)))
  categorical_data = data_modified$covid_status

  #finding number of locations
  num_loc = n_distinct(key)
  
  #length of data
  n = dim(regress_matrix)[1]
  
  num_time = n / num_loc
  #number of regressors
  len = dim(regress_matrix)[2]
  
  #chain initialisation for the three chains corresponding to Gelman-Rubin-Mechanism
  #chain1
  pi_sample = matrix(rep(0, n), nrow = n)
  omega_sample = matrix(rep(0, n), nrow = n)
  theta_sample = matrix(rep(0, len), nrow = len)
  sigma_w_sample = c(1)
  sigma_e_sample = c(1)
  
  pi0 = matrix(rep(0, n), nrow = n)[, 1]
  omega0 = matrix(rep(0, n), nrow = n)[, 1]
  theta0 = matrix(rep(0, len), nrow = len)[, 1]
  sigma_w0 = c(1)[1]
  sigma_e0 = c(1)[1]
  
  
  #chain2
  pi_sample1 = matrix(rep(0.4, n), nrow = n)
  omega_sample1 = matrix(rep(0.4, n), nrow = n)
  theta_sample1 = matrix(rep(0, len), nrow = len) + 0.1
  sigma_w_sample1 = c(0.7)
  sigma_e_sample1 = c(0.7)
  
  pi10 = matrix(rep(.4, n), nrow = n)[, 1]
  omega10 = matrix(rep(.4, n), nrow = n)[, 1]
  theta10 =  matrix(rep(0, len), nrow = len) + 0.1
  sigma_w10 = c(0.7)[1]
  sigma_e10 = c(0.7)[1]
  
  #chain3
  pi_sample2 = matrix(rep(0.9, n), nrow = n)
  omega_sample2 = matrix(rep(0.25, n), nrow = n)
  theta_sample2 = matrix(rep(0, len), nrow = len) + 0.2
  sigma_w_sample2 = c(1)
  sigma_e_sample2 = c(2)
  
  pi20 = matrix(rep(.9, n), nrow = n)[, 1]
  omega20 = matrix(rep(.25, n), nrow = n)[, 1]
  theta20 =  matrix(rep(0, len), nrow = len) + 0.2
  sigma_w20 = c(1)[1]
  sigma_e20 = c(2)[1]
  
  #Phi_initiate
  phi_sp1 = 1
  phi_sp2 = 2
  phi_sp3 = 3
  phi_tmp1 = 0.1
  phi_tmp2 = 0.2
  phi_tmp3 = 0.3
  
  i = 1
  
  value = 5
  
  count = 0
  
  while (value > GR_threshold & count < max_burn) {
    # start_time<-Sys.time()
    #########################################################################################33
    #Getting time difference matrix
    time_distance_mat = as.matrix(dist(time_points_common_used,
                                       diag  =  TRUE,
                                       upper  =  TRUE))
    tmp = exp(-phi_tmp * time_distance_mat)
    
    #Getting distance matrix
    distance_mat = as.matrix(dist(unique(cbind(
      location_mat[, c("lat")], location_mat[, c("long")]
    ))))###[seq(1,124,by=2),]###
    sp = exp(-phi_sp * distance_mat)
    
    #INTERCHANGE ALGORITHM
    interchange = list()
    for (j in 1:num_loc) {
      interchange[[j]] = c(1:num_loc)
      
      if (j != 1) {
        interchange[[j]][1] = j
        interchange[[j]][j] = 1
      }
    }
    
    #Getting inv sp_22 to use for every iteration
    
    mat_21 = vector()
    for (j in (1:num_loc)) {
      mat_21 = cbind(mat_21, (diag(num_loc)[interchange[[j]],] %*% sp %*% diag(num_loc)[interchange[[j]],])[-1, 1])
    }
    
    #counter for number of iterations that have occured
    count = 0
    
    inv_sp = chol2inv(chol(sp))
    inv_tmp = chol2inv(chol(tmp))
    
    
    
    
    #Getting inv_sp_22 to use for every iteration
    Inv_22 = list()
    for (j in (1:num_loc)) {
      inv_s_t = diag(num_loc)[interchange[[j]],] %*% inv_sp %*% diag(num_loc)[interchange[[j]],]
      inv_s_t_22 = inv_s_t[-1,-1] - (inv_s_t[-1, 1] %*% t(inv_s_t[1,-1])) /
        inv_s_t[1, 1]
      Inv_22[[j]] = inv_s_t_22
    }
    
    
    #following lists(a-d) have been introduced for ease in calculations in the gibbs sampling procedures
    
    #a
    sig = list()
    for (j in 1:num_loc) {
      sig[[j]] = (1 - t(mat_21[, j]) %*% Inv_22[[j]]
                  %*% mat_21[, j]) #%x% tmp
    }
    
    #b
    inv_sig1 = list()
    for (j in 1:num_loc) {
      inv_sig1[[j]] = chol2inv(chol(sig[[j]]))
    }
    
    #c
    inv_sig = list()
    for (j in 1:num_loc) {
      inv_sig[[j]] = chol2inv(chol((1 - t(mat_21[, j]) %*% Inv_22[[j]]
                                    %*% mat_21[, j]) %x% tmp))
    }
    
    #d
    mat_inv = list()
    for (j in 1:num_loc) {
      mat_inv[[j]] = ((t(mat_21[, j]) %*% Inv_22[[j]]))
    }
    
    #computed before hand for ease in calculations
    XtX_mat = solve(t(regress_matrix) %*% regress_matrix)
    
    XtXX_mat = XtX_mat %*% t(regress_matrix)
    ##########################################################################################################
    pi = pi0
    theta = theta0
    sigma_e = sigma_e0
    sigma_w = sigma_w0
    omega_mat = matrix(omega0, nrow = num_time)
    phi_sp = phi_sp1
    phi_tmp = phi_tmp1
    
    pi_mat = matrix(pi, nrow = num_time, ncol = num_loc)
    
    Regress_mat_theta = (regress_matrix %*% theta)
    Regress_mat_theta_mat = matrix(Regress_mat_theta, nrow = num_time, ncol = num_loc)
    
    #omega
    for (j in (1:num_loc)) {
      omega_interchange = matrix((omega_mat %*% diag(num_loc)[interchange[[j]],])
                                 [(1 + num_time):(n)] , nrow = num_time)
      
      mu_cond = c(diag(num_time) %*% omega_interchange %*% t(mat_inv[[j]]))
      
      inv = c(t(inv_tmp) %*% mu_cond %*% inv_sig1[[j]])
      
      inv_2 = chol2inv(chol(inv_sig[[j]] / sigma_w  + diag(num_time) / sigma_e))
      
      inv2_inv = inv_2 %*% (inv / sigma_w +
                              (pi_mat[, j] - Regress_mat_theta_mat[, j]) / sigma_e)
      
      omega_mat[, j] = inv2_inv + t(chol(inv_2))  %*% rnorm(num_time)
      
    }
    #theta
    theta_chol = t(chol(XtX_mat / sigma_e))
    theta =  XtXX_mat %*% (pi - omega_mat[1:n]) + theta_chol %*%
      rnorm(len)
    
    #sigma_w
    sigma_w_kro = sum(omega_mat[1:n] * c(inv_tmp %*% omega_mat %*% inv_sp))
    sigma_w = 1 / rgamma(1,
                         inv_gamma_parameter_1 + n / 2,
                         rate = sigma_w_kro / 2 + inv_gamma_parameter_2)
    
    #sigma_e
    D <- Regress_mat_theta + omega_mat[1:n]
    pi_minus_D_vector = pi - D
    pi_t_pi = t(pi_minus_D_vector) %*% (pi_minus_D_vector) / 2
    sigma_e = 1 / rgamma(1,
                         inv_gamma_parameter_1 + n / 2,
                         rate = (pi_t_pi + inv_gamma_parameter_2))
    
    
    #pi
    
    pi1 = sapply(
      D,
      FUN = function(x)
        rtruncnorm(1, 0, Inf, x, sqrt(sigma_e))
    )
    pi2 = sapply(
      D,
      FUN = function(x)
        rtruncnorm(1,-Inf, 0, x, sqrt(sigma_e))
    )
    u_mat_training = t(omega_mat)
    
    pi = pi1 * as.numeric(categorical_data == 1) + pi2 * as.numeric(categorical_data == 0)
    
    phi_v_sp_slice = diversitree::mcmc(
      lik = phi_sp_log_density_func,
      x.init = c(phi_sp),
      nsteps = 1,
      w = 5,
      lower = 0,
      upper = 3
    )
    
    phi_sp = phi_v_sp_slice$pars
    
    phi_v_tmp_slice = diversitree::mcmc(
      lik = phi_tmp_log_density_func,
      x.init = c(phi_tmp),
      nsteps = 1,
      w = 5,
      lower = 0,
      upper = 3
    )
    
    phi_tmp = phi_v_tmp_slice$pars
    
    
    
    
    {
      pi0 =  pi
      omega0 = omega_mat[1:n]
      theta0 = theta
      sigma_e0 = sigma_e
      sigma_w0 = sigma_w
      phi_sp1 = phi_sp
      phi_tmp1 = phi_tmp
    }
    
    #Collecting samples for prediction
    if (i %% gap == 0) {
      theta_sample = cbind(theta_sample, theta)
      sigma_e_sample = c(sigma_e_sample, sigma_e)
      sigma_w_sample = c(sigma_w_sample, sigma_w)
    }
    
    
    #CHAIN2
    #########################################################################################33
    #Getting time difference matrix
    time_distance_mat = as.matrix(dist(time_points_common_used,
                                       diag  =  TRUE,
                                       upper  =  TRUE))
    tmp = exp(-phi_tmp * time_distance_mat)
    
    #Getting distance matrix
    distance_mat = as.matrix(dist(unique(cbind(
      location_mat[, c("lat")], location_mat[, c("long")]
    ))))###[seq(1,124,by=2),]###
    sp = exp(-phi_sp * distance_mat)
    
    #INTERCHANGE ALGORITHM
    interchange = list()
    for (j in 1:num_loc) {
      interchange[[j]] = c(1:num_loc)
      
      if (j != 1) {
        interchange[[j]][1] = j
        interchange[[j]][j] = 1
      }
    }
    
    #Getting inv sp_22 to use for every iteration
    
    mat_21 = vector()
    for (j in (1:num_loc)) {
      mat_21 = cbind(mat_21, (diag(num_loc)[interchange[[j]],] %*% sp %*% diag(num_loc)[interchange[[j]],])[-1, 1])
    }
    
    #counter for number of iterations that have occured
    count = 0
    
    inv_sp = chol2inv(chol(sp))
    inv_tmp = chol2inv(chol(tmp))
    
    
    
    
    #Getting inv_sp_22 to use for every iteration
    Inv_22 = list()
    for (j in (1:num_loc)) {
      inv_s_t = diag(num_loc)[interchange[[j]],] %*% inv_sp %*% diag(num_loc)[interchange[[j]],]
      inv_s_t_22 = inv_s_t[-1,-1] - (inv_s_t[-1, 1] %*% t(inv_s_t[1,-1])) /
        inv_s_t[1, 1]
      Inv_22[[j]] = inv_s_t_22
    }
    
    
    #following lists(a-d) have been introduced for ease in calculations in the gibbs sampling procedures
    
    #a
    sig = list()
    for (j in 1:num_loc) {
      sig[[j]] = (1 - t(mat_21[, j]) %*% Inv_22[[j]]
                  %*% mat_21[, j]) #%x% tmp
    }
    
    #b
    inv_sig1 = list()
    for (j in 1:num_loc) {
      inv_sig1[[j]] = chol2inv(chol(sig[[j]]))
    }
    
    #c
    inv_sig = list()
    for (j in 1:num_loc) {
      inv_sig[[j]] = chol2inv(chol((1 - t(mat_21[, j]) %*% Inv_22[[j]]
                                    %*% mat_21[, j]) %x% tmp))
    }
    
    #d
    mat_inv = list()
    for (j in 1:num_loc) {
      mat_inv[[j]] = ((t(mat_21[, j]) %*% Inv_22[[j]]))
    }
    
    #computed before hand for ease in calculations
    XtX_mat = solve(t(regress_matrix) %*% regress_matrix)
    
    XtXX_mat = XtX_mat %*% t(regress_matrix)
    ##########################################################################################################
    #omega posterior
    pi = pi10
    theta = theta10
    sigma_e = sigma_e10
    sigma_w = sigma_w10
    omega_mat = matrix(omega10, nrow = num_time)
    phi_sp = phi_sp2
    phi_tmp = phi_tmp2
    
    pi_mat = matrix(pi, nrow = num_time, ncol = num_loc)
    
    Regress_mat_theta = (regress_matrix %*% theta)
    Regress_mat_theta_mat = matrix(Regress_mat_theta, nrow = num_time, ncol = num_loc)
    
    
    #omega
    for (j in (1:num_loc)) {
      omega_interchange = matrix((omega_mat %*% diag(num_loc)[interchange[[j]],])
                                 [(1 + num_time):(n)] , nrow = num_time)
      
      mu_cond = c(diag(num_time) %*% omega_interchange %*% t(mat_inv[[j]]))
      
      inv = c(t(inv_tmp) %*% mu_cond %*% inv_sig1[[j]])
      
      inv_2 = chol2inv(chol(inv_sig[[j]] / sigma_w  + diag(num_time) / sigma_e))
      
      inv2_inv = inv_2 %*% (inv / sigma_w +
                              (pi_mat[, j] - Regress_mat_theta_mat[, j]) / sigma_e)
      
      omega_mat[, j] = inv2_inv + t(chol(inv_2))  %*% rnorm(num_time)
      
    }
    #theta
    theta_chol = t(chol(XtX_mat / sigma_e))
    theta =  XtXX_mat %*% (pi - omega_mat[1:n]) + theta_chol %*%
      rnorm(len)
    
    #sigma_w
    sigma_w_kro = sum(omega_mat[1:n] * c(inv_tmp %*% omega_mat %*% inv_sp))
    sigma_w = 1 / rgamma(1,
                         inv_gamma_parameter_1 + n / 2,
                         rate = sigma_w_kro / 2 + inv_gamma_parameter_2)
    
    #sigma_e
    D <- Regress_mat_theta + omega_mat[1:n]
    pi_minus_D_vector = pi - D
    pi_t_pi = t(pi_minus_D_vector) %*% (pi_minus_D_vector) / 2
    sigma_e = 1 / rgamma(1,
                         inv_gamma_parameter_1 + n / 2,
                         rate = (pi_t_pi + inv_gamma_parameter_2))
    
    
    #pi
    
    pi1 = sapply(
      D,
      FUN = function(x)
        rtruncnorm(1, 0, Inf, x, sqrt(sigma_e))
    )
    pi2 = sapply(
      D,
      FUN = function(x)
        rtruncnorm(1,-Inf, 0, x, sqrt(sigma_e))
    )
    u_mat_training = t(omega_mat)
    
    pi = pi1 * as.numeric(categorical_data == 1) + pi2 * as.numeric(categorical_data == 0)
    
    phi_v_sp_slice = diversitree::mcmc(
      lik = phi_sp_log_density_func,
      x.init = c(phi_sp),
      nsteps = 1,
      w = 5,
      lower = 0,
      upper = 3
    )
    
    phi_sp = phi_v_sp_slice$pars
    
    phi_v_tmp_slice = diversitree::mcmc(
      lik = phi_tmp_log_density_func,
      x.init = c(phi_tmp),
      nsteps = 1,
      w = 5,
      lower = 0,
      upper = 3
    )
    
    phi_tmp = phi_v_tmp_slice$pars
    
    
    
    {
      pi10 =  pi
      omega10 = omega_mat[1:n]
      theta10 = theta
      sigma_e10 = sigma_e
      sigma_w10 = sigma_w
      phi_sp2 = phi_sp
      phi_tmp2 = phi_tmp
    }
    
    
    #Collecting samples for prediction
    if (i %% gap == 0) {
      theta_sample1 = cbind(theta_sample1, theta)
      sigma_e_sample1 = c(sigma_e_sample1, sigma_e)
      
      sigma_w_sample1 = c(sigma_w_sample1, sigma_w)
    }
    
    
    
    #CHAIN3
    #########################################################################################33
    #Getting time difference matrix
    time_distance_mat = as.matrix(dist(time_points_common_used,
                                       diag  =  TRUE,
                                       upper  =  TRUE))
    tmp = exp(-phi_tmp * time_distance_mat)
    
    #Getting distance matrix
    distance_mat = as.matrix(dist(unique(cbind(
      location_mat[, c("lat")], location_mat[, c("long")]
    ))))###[seq(1,124,by=2),]###
    sp = exp(-phi_sp * distance_mat)
    
    #INTERCHANGE ALGORITHM
    interchange = list()
    for (j in 1:num_loc) {
      interchange[[j]] = c(1:num_loc)
      
      if (j != 1) {
        interchange[[j]][1] = j
        interchange[[j]][j] = 1
      }
    }
    
    #Getting inv sp_22 to use for every iteration
    
    mat_21 = vector()
    for (j in (1:num_loc)) {
      mat_21 = cbind(mat_21, (diag(num_loc)[interchange[[j]],] %*% sp %*% diag(num_loc)[interchange[[j]],])[-1, 1])
    }
    
    #counter for number of iterations that have occured
    count = 0
    
    inv_sp = chol2inv(chol(sp))
    inv_tmp = chol2inv(chol(tmp))
    
    
    
    
    #Getting inv_sp_22 to use for every iteration
    Inv_22 = list()
    for (j in (1:num_loc)) {
      inv_s_t = diag(num_loc)[interchange[[j]],] %*% inv_sp %*% diag(num_loc)[interchange[[j]],]
      inv_s_t_22 = inv_s_t[-1,-1] - (inv_s_t[-1, 1] %*% t(inv_s_t[1,-1])) /
        inv_s_t[1, 1]
      Inv_22[[j]] = inv_s_t_22
    }
    
    
    #following lists(a-d) have been introduced for ease in calculations in the gibbs sampling procedures
    
    #a
    sig = list()
    for (j in 1:num_loc) {
      sig[[j]] = (1 - t(mat_21[, j]) %*% Inv_22[[j]]
                  %*% mat_21[, j]) #%x% tmp
    }
    
    #b
    inv_sig1 = list()
    for (j in 1:num_loc) {
      inv_sig1[[j]] = chol2inv(chol(sig[[j]]))
    }
    
    #c
    inv_sig = list()
    for (j in 1:num_loc) {
      inv_sig[[j]] = chol2inv(chol((1 - t(mat_21[, j]) %*% Inv_22[[j]]
                                    %*% mat_21[, j]) %x% tmp))
    }
    
    #d
    mat_inv = list()
    for (j in 1:num_loc) {
      mat_inv[[j]] = ((t(mat_21[, j]) %*% Inv_22[[j]]))
    }
    
    #computed before hand for ease in calculations
    XtX_mat = solve(t(regress_matrix) %*% regress_matrix)
    
    XtXX_mat = XtX_mat %*% t(regress_matrix)
    ##########################################################################################################
    pi = pi20
    theta = theta20
    sigma_e = sigma_e20
    sigma_w = sigma_w20
    omega_mat = matrix(omega20, nrow = num_time)
    phi_sp = phi_sp3
    phi_tmp = phi_tmp3
    
    
    pi_mat = matrix(pi, nrow = num_time, ncol = num_loc)
    
    Regress_mat_theta = (regress_matrix %*% theta)
    Regress_mat_theta_mat = matrix(Regress_mat_theta, nrow = num_time, ncol = num_loc)
    
    #omega
    for (j in (1:num_loc)) {
      omega_interchange = matrix((omega_mat %*% diag(num_loc)[interchange[[j]],])
                                 [(1 + num_time):(n)] , nrow = num_time)
      
      mu_cond = c(diag(num_time) %*% omega_interchange %*% t(mat_inv[[j]]))
      
      inv = c(t(inv_tmp) %*% mu_cond %*% inv_sig1[[j]])
      
      inv_2 = chol2inv(chol(inv_sig[[j]] / sigma_w  + diag(num_time) / sigma_e))
      
      inv2_inv = inv_2 %*% (inv / sigma_w +
                              (pi_mat[, j] - Regress_mat_theta_mat[, j]) / sigma_e)
      
      omega_mat[, j] = inv2_inv + t(chol(inv_2))  %*% rnorm(num_time)
      
    }
    #theta
    theta_chol = t(chol(XtX_mat / sigma_e))
    theta =  XtXX_mat %*% (pi - omega_mat[1:n]) + theta_chol %*%
      rnorm(len)
    
    #sigma_w
    sigma_w_kro = sum(omega_mat[1:n] * c(inv_tmp %*% omega_mat %*% inv_sp))
    sigma_w = 1 / rgamma(1,
                         inv_gamma_parameter_1 + n / 2,
                         rate = sigma_w_kro / 2 + inv_gamma_parameter_2)
    
    #sigma_e
    D <- Regress_mat_theta + omega_mat[1:n]
    pi_minus_D_vector = pi - D
    pi_t_pi = t(pi_minus_D_vector) %*% (pi_minus_D_vector) / 2
    sigma_e = 1 / rgamma(1,
                         inv_gamma_parameter_1 + n / 2,
                         rate = (pi_t_pi + inv_gamma_parameter_2))
    
    
    #pi
    
    pi1 = sapply(
      D,
      FUN = function(x)
        rtruncnorm(1, 0, Inf, x, sqrt(sigma_e))
    )
    pi2 = sapply(
      D,
      FUN = function(x)
        rtruncnorm(1,-Inf, 0, x, sqrt(sigma_e))
    )
    u_mat_training = t(omega_mat)
    
    pi = pi1 * as.numeric(categorical_data == 1) + pi2 * as.numeric(categorical_data == 0)
    
    phi_v_sp_slice = diversitree::mcmc(
      lik = phi_sp_log_density_func,
      x.init = c(phi_sp),
      nsteps = 1,
      w = 5,
      lower = 0,
      upper = 3
    )
    
    phi_sp = phi_v_sp_slice$pars
    
    phi_v_tmp_slice = diversitree::mcmc(
      lik = phi_tmp_log_density_func,
      x.init = c(phi_tmp),
      nsteps = 1,
      w = 5,
      lower = 0,
      upper = 3
    )
    
    phi_tmp = phi_v_tmp_slice$pars
    
    
    
    {
      pi20 =  pi
      omega20 = omega_mat[1:n]
      theta20 = theta
      sigma_e20 = sigma_e
      sigma_w20 = sigma_w
      phi_sp3 = phi_sp
      phi_tmp3 = phi_tmp
    }
    
    
    #Collecting samples for prediction
    if (i %% gap == 0) {
      theta_sample2 = cbind(theta_sample2, theta)
      sigma_e_sample2 = c(sigma_e_sample2, sigma_e)
      
      sigma_w_sample2 = c(sigma_w_sample2, sigma_w)
    }
    
    
    count = count + 1
    print(count)
    
    if (i > 100 && i %% gap == 0) {
      p_mcmc_obj = list()
      for (j in 1:(i / gap)) {
        p_mcmc_obj[[j]] = ((cbind(
          theta_sample[, j],
          theta_sample1[, j],
          theta_sample2[, j]
        )))
      }
  
      plist = lapply(apply(do.call(
        rbind, lapply(p_mcmc_obj, function(x)
          x[1,])
      ), 2, as.data.frame), coda::mcmc)
      
      
      value = as.numeric((gelman.diag(x = plist)$psrf)[1])
      
    }
  }
  
  
 
  
  print("the code has converged YAAAYYYY")
  theta_sample4 = c()
  sigma_e_sample4 = c()
  sigma_w_sample4 = c()
  pi_sample4 = c()
  omega_sample4 = c()
  phi_sp_store = c()
  phi_tmp_store = c()
  
  for (l in 1:data_modified_size) {
    ##########################################################################################33
    #Getting time difference matrix
    time_distance_mat = as.matrix(dist(time_points_common_used,
                                       diag  =  TRUE,
                                       upper  =  TRUE))
    tmp = exp(-phi_tmp * time_distance_mat)
    
    #Getting distance matrix
    distance_mat = as.matrix(dist(unique(cbind(
      location_mat[, c("lat")], location_mat[, c("long")]
    ))))###[seq(1,124,by=2),]###
    sp = exp(-phi_sp * distance_mat)
    
    #INTERCHANGE ALGORITHM
    interchange = list()
    for (j in 1:num_loc) {
      interchange[[j]] = c(1:num_loc)
      
      if (j != 1) {
        interchange[[j]][1] = j
        interchange[[j]][j] = 1
      }
    }
    
    #Getting inv sp_22 to use for every iteration
    
    mat_21 = vector()
    for (j in (1:num_loc)) {
      mat_21 = cbind(mat_21, (diag(num_loc)[interchange[[j]],] %*% sp %*% diag(num_loc)[interchange[[j]],])[-1, 1])
    }
    
    #counter for number of iterations that have occured
    count = 0
    
    inv_sp = chol2inv(chol(sp))
    inv_tmp = chol2inv(chol(tmp))
    
    
    
    
    #Getting inv_sp_22 to use for every iteration
    Inv_22 = list()
    for (j in (1:num_loc)) {
      inv_s_t = diag(num_loc)[interchange[[j]],] %*% inv_sp %*% diag(num_loc)[interchange[[j]],]
      inv_s_t_22 = inv_s_t[-1,-1] - (inv_s_t[-1, 1] %*% t(inv_s_t[1,-1])) /
        inv_s_t[1, 1]
      Inv_22[[j]] = inv_s_t_22
    }
    
    
    #following lists(a-d) have been introduced for ease in calculations in the gibbs sampling procedures
    
    #a
    sig = list()
    for (j in 1:num_loc) {
      sig[[j]] = (1 - t(mat_21[, j]) %*% Inv_22[[j]]
                  %*% mat_21[, j]) #%x% tmp
    }
    
    #b
    inv_sig1 = list()
    for (j in 1:num_loc) {
      inv_sig1[[j]] = chol2inv(chol(sig[[j]]))
    }
    
    #c
    inv_sig = list()
    for (j in 1:num_loc) {
      inv_sig[[j]] = chol2inv(chol((1 - t(mat_21[, j]) %*% Inv_22[[j]]
                                    %*% mat_21[, j]) %x% tmp))
    }
    
    #d
    mat_inv = list()
    for (j in 1:num_loc) {
      mat_inv[[j]] = ((t(mat_21[, j]) %*% Inv_22[[j]]))
    }
    
    #computed before hand for ease in calculations
    XtX_mat = solve(t(regress_matrix) %*% regress_matrix)
    
    XtXX_mat = XtX_mat %*% t(regress_matrix)
    ###########################################################################################################
 
    #omega posterior
    
    pi_mat = matrix(pi, nrow = num_time, ncol = num_loc)
    
    Regress_mat_theta = (regress_matrix %*% theta)
    Regress_mat_theta_mat = matrix(Regress_mat_theta, nrow = num_time, ncol = num_loc)
    
    #omega
    for (j in (1:num_loc)) {
      omega_interchange = matrix((omega_mat %*% diag(num_loc)[interchange[[j]],])
                                 [(1 + num_time):(n)] , nrow = num_time)
      
      mu_cond = c(diag(num_time) %*% omega_interchange %*% t(mat_inv[[j]]))
      
      inv = c(t(inv_tmp) %*% mu_cond %*% inv_sig1[[j]])
      
      inv_2 = chol2inv(chol(inv_sig[[j]] / sigma_w  + diag(num_time) / sigma_e))
      
      inv2_inv = inv_2 %*% (inv / sigma_w +
                              (pi_mat[, j] - Regress_mat_theta_mat[, j]) / sigma_e)
      
      omega_mat[, j] = inv2_inv + t(chol(inv_2))  %*% rnorm(num_time)
      
    }
    #theta
    theta_chol = t(chol(XtX_mat / sigma_e))
    theta =  XtXX_mat %*% (pi - omega_mat[1:n]) + theta_chol %*%
      rnorm(len)
    
    #sigma_w
    sigma_w_kro = sum(omega_mat[1:n] * c(inv_tmp %*% omega_mat %*% inv_sp))
    sigma_w = 1 / rgamma(1,
                         inv_gamma_parameter_1 + n / 2,
                         rate = sigma_w_kro / 2 + inv_gamma_parameter_2)
    
    #sigma_e
    D <- Regress_mat_theta + omega_mat[1:n]
    pi_minus_D_vector = pi - D
    pi_t_pi = t(pi_minus_D_vector) %*% (pi_minus_D_vector) / 2
    sigma_e = 1 / rgamma(1,
                         inv_gamma_parameter_1 + n / 2,
                         rate = (pi_t_pi + inv_gamma_parameter_2))
    
    
    #pi
    
    pi1 = sapply(
      D,
      FUN = function(x)
        rtruncnorm(1, 0, Inf, x, sqrt(sigma_e))
    )
    pi2 = sapply(
      D,
      FUN = function(x)
        rtruncnorm(1,-Inf, 0, x, sqrt(sigma_e))
    )
    
    u_mat_training = t(omega_mat)
    
    pi = pi1 * as.numeric(categorical_data == 1) + pi2 * as.numeric(categorical_data == 0)
    
    phi_v_sp_slice = diversitree::mcmc(
      lik = phi_sp_log_density_func,
      x.init = c(phi_sp),
      nsteps = 1,
      w = 5,
      lower = 0,
      upper = 3
    )
    
    phi_sp = phi_v_sp_slice$pars
    
    phi_v_tmp_slice = diversitree::mcmc(
      lik = phi_tmp_log_density_func,
      x.init = c(phi_tmp),
      nsteps = 1,
      w = 5,
      lower = 0,
      upper = 3
    )
    
    phi_tmp = phi_v_tmp_slice$pars
    
    print(theta)
    
    
    #Collection of samples
    
    if (l %% gap == 0) {
      #pi_sample4 = cbind(pi_sample4, pi)
      omega_sample4 = cbind(omega_sample4, omega_mat[1:n])
      theta_sample4 = cbind(theta_sample4, theta)
      sigma_e_sample4 = c(sigma_e_sample4, sigma_e)
      sigma_w_sample4 = c(sigma_w_sample4, sigma_w)
      phi_sp_store = c(phi_sp_store, phi_sp)
      phi_tmp_store = c(phi_tmp_store, phi_tmp)
      # PLOT = plot.ts(t(rbind(
      #   #sigma_w_sample4,
      #   sigma_e_sample4,theta_sample4[c(1:len),]
      # )))
    }
    
    
    print(count)
    count = count + 1
    l = l + 1
  }
