# Load required packages
library(tidyverse)
library(dplyr)
library(truncnorm)
library(coda)
library(pROC)
library(ROCR)
library(measures)
library(geosphere)
library(ProbitSpatial)
library(bayestestR)
library(sp)
library(spdep)
library(spatialreg)
library(CARBayes)
library(e1071)

# Parameters
phi_sp <- 1
phi_tmp <- 0.1
data_modified_size <- 1000
gap <- 10
GR_threshold <- 1.1
max_burn <- 10000
inv_gamma_parameter_1 <- 2
inv_gamma_parameter_2 <- 1
l1 <- 0
l2 <- 0

# Regression formula
reg_f <- "covid_status ~ prop_vaccinated + I(log(death+1)) + lagged_death + I(log(population)) +
  linear_trend + quadratic_trend + sine_trend + cosine_trend"
mmat <- model.matrix(formula(reg_f), data = dff)
regress_matrix <- mmat

# Categorical data
categorical_data <- dff$covid_status 

# Length of data
n <- dim(regress_matrix)[1]

# Number of locations and time points
num_loc <- n_distinct(key)
num_time <- n / num_loc

# Number of regressors
len <- dim(regress_matrix)[2]

# Initial values from GLM fit
glm_train <- glm(
  formula = formula(reg_f),
  family = binomial(link = "probit"),
  data = dff
)
initvalue <- glm_train$coefficients

# Precomputed matrices
XtX_mat <- solve(t(regress_matrix) %*% regress_matrix)
XtXX_mat <- XtX_mat %*% t(regress_matrix)

# Iteration loop
for (l in 1:data_modified_size) {
  
  # Compute exponential of distance matrices
  tmp <- exp(-phi_tmp * time_distance_mat)
  sp <- exp(-phi_sp * distance_mat)
  
  # Interchange algorithm
  interchange <- list()
  for (j in 1:num_loc) {
    interchange[[j]] <- c(1:num_loc)
    if (j != 1) {
      interchange[[j]][1] <- j
      interchange[[j]][j] <- 1
    }
  }
  
  # Compute inv_sp_22
  mat_21 <- vector()
  for (j in 1:num_loc) {
    mat_21 <- cbind(
      mat_21,
      (diag(num_loc)[interchange[[j]], ] %*% sp %*% diag(num_loc)[interchange[[j]], ])[-1, 1]
    )
  }
  inv_sp <- chol2inv(chol(sp))
  inv_tmp <- chol2inv(chol(tmp))
  
  # Compute Inv_22
  Inv_22 <- list()
  for (j in 1:num_loc) {
    inv_s_t <- diag(num_loc)[interchange[[j]], ] %*% inv_sp %*% diag(num_loc)[interchange[[j]], ]
    inv_s_t_22 <- inv_s_t[-1, -1] - (inv_s_t[-1, 1] %*% t(inv_s_t[1, -1])) / inv_s_t[1, 1]
    Inv_22[[j]] <- inv_s_t_22
  }
  
  # Lists for ease in calculations
  # a
  sig <- list()
  for (j in 1:num_loc) {
    sig[[j]] <- (1 - t(mat_21[, j]) %*% Inv_22[[j]] %*% mat_21[, j])
  }
  
  # b
  inv_sig1 <- list()
  for (j in 1:num_loc) {
    inv_sig1[[j]] <- chol2inv(chol(sig[[j]]))
  }
  
  # c
  inv_sig <- list()
  for (j in 1:num_loc) {
    inv_sig[[j]] <- chol2inv(chol((1 - t(mat_21[, j]) %*% Inv_22[[j]] %*% mat_21[, j]) %x% tmp))
  }
  
  # d
  mat_inv <- list()
  for (j in 1:num_loc) {
    mat_inv[[j]] <- (t(mat_21[, j]) %*% Inv_22[[j]])
  }
  
  # Initialize matrices
  pi_mat <- matrix(pi, nrow = num_time, ncol = num_loc)
  Regress_mat_theta <- (regress_matrix %*% theta)
  Regress_mat_theta_mat <- matrix(Regress_mat_theta, nrow = num_time, ncol = num_loc)
  
  # Omega
  for (j in 1:num_loc) {
    omega_interchange <- matrix((omega_mat %*% diag(num_loc)[interchange[[j]], ])
                                [(1 + num_time):(n)], nrow = num_time)
    mu_cond <- c(diag(num_time) %*% omega_interchange %*% t(mat_inv[[j]]))
    inv <- c(t(inv_tmp) %*% mu_cond %*% inv_sig1[[j]])
    inv_2 <- chol2inv(chol(inv_sig[[j]] / sigma_w  + diag(num_time) / sigma_e))
    inv2_inv <- inv_2 %*% (inv / sigma_w + (pi_mat[, j] - Regress_mat_theta_mat[, j]) / sigma_e)
    omega_mat[, j] <- inv2_inv + t(chol(inv_2)) %*% rnorm(num_time)
  }
  
  # Theta
  theta_chol <- t(chol(XtX_mat / sigma_e))
  theta <-  XtXX_mat %*% (pi - omega_mat[1:n]) + theta_chol %*% rnorm(len)
  
  # sigma_w
  sigma_w_kro <- sum(omega_mat[1:n] * c(inv_tmp %*% omega_mat %*% inv_sp))
  sigma_w <- 1 / rgamma(1,
                        inv_gamma_parameter_1 + n / 2,
                        rate = sigma_w_kro / 2 + inv_gamma_parameter_2)
  
  # sigma_e
  D <- Regress_mat_theta + omega_mat[1:n]
  pi_minus_D_vector <- pi - D
  pi_t_pi <- t(pi_minus_D_vector) %*% (pi_minus_D_vector) / 2
  sigma_e <- 1 / rgamma(1,
                        inv_gamma_parameter_1 + n / 2,
                        rate = (pi_t_pi + inv_gamma_parameter_2))
  
  # pi
  pi1 <- sapply(
    D,
    FUN = function(x)
      rtruncnorm(1, 0, Inf, x, sqrt(sigma_e))
  )
  pi2 <- sapply(
    D,
    FUN = function(x)
      rtruncnorm(1, -Inf, 0, x, sqrt(sigma_e))
  )
  
  u_mat_training <- t(omega_mat)
  pi <- pi1 * as.numeric(categorical_data == 1) + pi2 * as.numeric(categorical_data == 0)
  
  # Phi_sp log density function
  phi_sp_log_density_func <- function(phi_sp) {
    SIGMA_sp <- exp(-phi_sp * distance_mat)
    log_det_SIGMA_sp <- (determinant(SIGMA_sp, logarithm = TRUE))$modulus[1]
    inv_SIGMA_sp <- chol2inv(chol(SIGMA_sp))
    mat_mult_inv_kro_val_training <- 0
    for (col1 in 1:num_time) {
      for (col2 in col1:num_time) {
        if (col1 == col2) {
          mat_mult_inv_kro_val_training <- (mat_mult_inv_kro_val_training +
                                              inv_tmp[col1, col2] *
                                              as.numeric(
                                                t(u_mat_training[, col2] %*%
                                                    inv_SIGMA_sp %*% u_mat_training[, col1])
                                              ))
        } else{
          mat_mult_inv_kro_val_training <- (mat_mult_inv_kro_val_training +
                                              (inv_tmp[col1, col2] + inv_tmp[col2, col1]) *
                                              as.numeric(
                                                t(u_mat_training[, col2] %*%
                                                    inv_SIGMA_sp %*% u_mat_training[, col1])
                                              ))
        }
      }
    }
    (-(as.numeric(mat_mult_inv_kro_val_training)) / (2 * sigma_w) -
        num_time * log_det_SIGMA_sp / 2) + dgamma(phi_sp, shape = 2, scale = 5, log = TRUE)
  }
  
  # Phi_sp MCMC sampling
  phi_v_sp_slice <- diversitree::mcmc(
    lik = phi_sp_log_density_func,
    x.init = c(phi_sp),
    nsteps = 1,
    w = 5,
    lower = l1,
    upper = 3,
    print.every = 0
  )
  
  phi_sp <- phi_v_sp_slice$pars
  
  # Phi_tmp log density function
  phi_tmp_log_density_func <- function(phi_tmp) {
    SIGMA_tmp <- exp(-phi_tmp * time_distance_mat)
    log_det_SIGMA_tmp <- (determinant(SIGMA_tmp, logarithm = TRUE))$modulus[1]
    inv_SIGMA_tmp <- chol2inv(chol(SIGMA_tmp))
    mat_mult_inv_kro_val2_training <- 0
    for (col1 in 1:num_time) {
      for (col2 in col1:num_time) {
        if (col1 == col2) {
          mat_mult_inv_kro_val2_training <- (mat_mult_inv_kro_val2_training +
                                               inv_SIGMA_tmp[col1, col2] *
                                               as.numeric(t(
                                                 u_mat_training[, col2] %*%
                                                   inv_sp %*% u_mat_training[, col1]
                                               )))
        } else{
          mat_mult_inv_kro_val2_training <- (
            mat_mult_inv_kro_val2_training +
              (inv_SIGMA_tmp[col1, col2] + inv_SIGMA_tmp[col2, col1]) *
              as.numeric(t(
                u_mat_training[, col2] %*%
                  inv_sp %*% u_mat_training[, col1]
              ))
          )
        }
      }
    }
    (-(as.numeric(mat_mult_inv_kro_val2_training)) / (2 * sigma_w) -
        num_loc * log_det_SIGMA_tmp / 2) + dgamma(phi_tmp, shape = 2, scale = 5, log = TRUE)
  }
  
  # Phi_tmp MCMC sampling
  phi_v_tmp_slice <- diversitree::mcmc(
    lik = phi_tmp_log_density_func,
    x.init = c(phi_tmp),
    nsteps = 1,
    w = 5,
    lower = l2,
    upper = 3,
    print.every = 0
  )
  
  phi_tmp <- phi_v_tmp_slice$pars
  
  # Collection of samples
  if (l %% gap == 0) {
    omega_sample4 <- cbind(omega_sample4, omega_mat[1:n])
    theta_sample4 <- cbind(theta_sample4, theta)
    sigma_e_sample4 <- c(sigma_e_sample4, sigma_e)
    sigma_w_sample4 <- c(sigma_w_sample4, sigma_w)
    phi_sp_store <- c(phi_sp_store, phi_sp)
    phi_tmp_store <- c(phi_tmp_store, phi_tmp)
  }
}
