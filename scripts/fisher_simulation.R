# ----------------------------------------------------
# Licensed under the Apache License 2.0.
# Written by Aakash Agrawal (akash.gagrawal@gmail.com)
# ----------------------------------------------------
# %
# Simulation for MOBW using FIM #

rm(list=ls())
options(warn=-1)

source("scripts/utils.R")
library(VGAM)

n_dpts = 400 # sample size
censoring = FALSE # whether to have censoring or not
tau = 0.5 # censoring factor

n_iter = 1000 # Number of monte carlo replications
true_param = c(2.35,1.11,1.92,1.63)

# Estimates
alpha_est = lam_1_est = lam_2_est = lam_3_est = MLE_est = rep(0,n_iter)

# CIs
ci_aph_l = ci_aph_u = ci_lam1_l = ci_lam1_u = rep(0,n_iter)
ci_lam2_l = ci_lam2_u = ci_lam3_l = ci_lam3_u = rep(0,n_iter)

# Average Length OF CIs.
AL_aph = AL_lam1 = AL_lam2 = AL_lam3 = 0 

# Coverage Percentage of CIs.
cp_aph = cp_lam1 = cp_lam2 = cp_lam3 = 0

# Convergence of MLEs
convergence = rep(0, n_iter)

# functions to get the probability factors
f_integral_weibull_1 <- function(t) {
  return(alpha * lam_1 * (t ** (lam_1 - 1)) * 
           exp(-alpha * ((t ** lam_1) + (t ** lam_2) + (t ** lam_3))))
}

f_integral_weibull_2 <- function(t) {
  return(alpha * lam_2 * (t ** (lam_2 - 1)) * 
           exp(-alpha * ((t ** lam_1) + (t ** lam_2) + (t ** lam_3))))
}

# Log Likelihood function
LL <- function(x){
  
  ## alpha, lambda_1, lambda_2, lambda_3 are the estimates we want to find
  sum_ = 0
  sum_j = 0
  
  alpha <<- x[1]
  lam_1 <<- x[2]
  lam_2 <<- x[3]
  lam_3 <<- x[4]
  
  p1 = integrate(f_integral_weibull_1, lower = 0, upper = Inf)$value
  p2 = integrate(f_integral_weibull_2, lower = 0, upper = Inf)$value
  p_star = (1 - p1 - p2)
  
  for(i in 1:nrow(lifetime_df)) {
    row <- lifetime_df[i, ]
    t = row$lifetime
    
    lambdas = c(lam_1, lam_2, lam_3)
    for (j in lambdas) {
      sum_j = sum_j + log((1 - G(t, j)))
    }
    
    # failure due to the first cause or treatment
    if (row$cause == 1) {
      sum_ = sum_ + log(-R_dash(t, lam_1, 0) / R(t, lam_1, 0))
    }
    
    # failure due to the second cause or lack of treatment
    if (row$cause == 2) {
      sum_ = sum_ + log(-R_dash(t, lam_2, 0) / R(t, lam_2, 0))
    }
    
    # failure due to both causes
    if (row$cause == 0) {
      sum_ = sum_ + 
        log(-(R_dash(t, lam_3, 0) / R(t, lam_3, 0)) - 
              (R_dash(t, lam_1, lam_2) / R(t, lam_1, lam_2)))
    }
  }
  
  sum_ = sum_ + alpha * sum_j + 
    (m_0 + m_1 + m_2) * log(alpha) + m_0 * log(p_star)
  
  return(sum_)
}

# Simulation Dataset
get_dataset <- function(x){
  # initialise a dataframe
  cols = c("cause", "lifetime")
  df = data.frame(matrix(nrow = n_dpts, ncol = length(cols)))
  
  # assign column names
  colnames(df) = cols

  u1 <- rweibull(n_dpts, shape = x[2], scale = 1 / (x[1] ** (1 / x[2])))
  u2 <- rweibull(n_dpts, shape = x[3], scale = 1 / (x[1] ** (1 / x[3])))
  u3 <- rweibull(n_dpts, shape = x[4], scale = 1 / (x[1] ** (1 / x[4])))
  u_bind = cbind(u1, u2, u3)
  
  df$cause = apply(u_bind, 1, which.min)
  df[df$cause == 3, ]$cause = 0 # tie case
  
  df$lifetime = apply(u_bind, 1, min)
  
  ## add some cases of censoring
  if (censoring) {
    df[df$lifetime > tau, ]$cause = 3
    df[df$lifetime > tau, ]$lifetime = tau
  }
  
  # set num fails in each cause
  set_num_fails(df)
  
  return(df)
}

i = 1
mle_error = 0
hessian_error = 0
start.time <- Sys.time()
while (i <= n_iter) {
  # track progress
  if(i %% 1 == 0){
    cat("iteration num: ", i)
    cat("\n")
  }
  
  # dataset
  df = get_dataset(true_param)
  print(c(m_0, m_1, m_2, m_c))
  
  lifetime_df <<- df # make df global
  
  # MLE
  MLE = try(
    {constrOptim(true_param,
                  LL,
                  grad = NULL,
                  control = list(fnscale=-1),
                  ui = diag(4),
                  ci = rep(0,4))},
    silent = TRUE
  )
  
  # if error then skip
  if(inherits(MLE, 'try-error')) {
    mle_error = mle_error + 1 
    next
  }

  # log the values
  MLE_est[i] = MLE$value
  alpha_est[i] = aph = MLE$par[1]
  lam_1_est[i] = l1 = MLE$par[2]
  lam_2_est[i] = l2 = MLE$par[3]
  lam_3_est[i] = l3 = MLE$par[4]
  
  cat("alpha_0", MLE$par[4] - true_param[4])
  cat("\n")
  
  # convergence 0 indicates successful completion 
  # 1 indicates that the iteration limit maxit had been reached.
  convergence[i] = conv = MLE$convergence
  
  # covariance matrix is inverse of negative of hessian
  CV = inv(-hessian(LL, c(aph, l1, l2, l3)))
  
  # if diagonal element less than 0
  if (any(diag(CV) < 0) | all(is.na(CV))) {
    hessian_error = hessian_error + 1
    next
  }
  
  # Divergence Cases
  if (abs( MLE$par[4] - true_param[4]) > 10) {
    next
  }
  
  # finding the 95% CI
  ci_aph_l[i] = aph - 1.96 * sqrt(CV[1, 1])
  ci_aph_u[i] = aph + 1.96 * sqrt(CV[1, 1])
  ci_lam1_l[i] = l1 - 1.96 * sqrt(CV[2, 2])
  ci_lam1_u[i] = l1 + 1.96 * sqrt(CV[2, 2])
  ci_lam2_l[i] = l2 - 1.96 * sqrt(CV[3, 3])
  ci_lam2_u[i] = l2 + 1.96 * sqrt(CV[3, 3])
  ci_lam3_l[i] = l3 - 1.96 * sqrt(CV[4, 4])
  ci_lam3_u[i] = l3 + 1.96 * sqrt(CV[4, 4])
  
  i = i + 1
}

end.time <- Sys.time()
time.taken <- end.time - start.time

cat("time taken for simulation: ", time.taken)

## Average Estimate
AE_aph = sum(alpha_est)/n_iter
AE_lam1 = sum(lam_1_est)/n_iter  
AE_lam2 = sum(lam_2_est)/n_iter
AE_lam3 = sum(lam_3_est)/n_iter

# MSE for each parameter
MSE_aph = sum((alpha_est - true_param[1])**2)/n_iter
MSE_lam1 = sum((lam_1_est - true_param[2])**2)/n_iter
MSE_lam2 = sum((lam_2_est - true_param[3])**2)/n_iter
MSE_lam3 = sum((lam_3_est - true_param[4])**2)/n_iter

# BIAS
BIAS_aph = AE_aph - true_param[1]
BIAS_lam1 = AE_lam1 - true_param[2]
BIAS_lam2 = AE_lam2 - true_param[3]
BIAS_lam3 = AE_lam3 - true_param[4]

# CP and AL
for (i in 1:n_iter) {
  AL_aph = AL_aph + ci_aph_u[i] - ci_aph_l[i]
  AL_lam1 = AL_lam1 + ci_lam1_u[i] - ci_lam1_l[i]
  AL_lam2 = AL_lam2 + ci_lam2_u[i] - ci_lam2_l[i]
  AL_lam3 = AL_lam3 + ci_lam3_u[i] - ci_lam3_l[i]
  
  if ((ci_aph_l[i] <= true_param[1]) && (ci_aph_u[i] >= true_param[1])) {
    cp_aph = cp_aph + 1
  }
  
  if((ci_lam1_l[i] <= true_param[2]) && (ci_lam1_u[i] >= true_param[2])){
    cp_lam1 = cp_lam1 + 1
  }
  
  if((ci_lam2_l[i] <= true_param[3]) && (ci_lam2_u[i] >= true_param[3])){
    cp_lam2 = cp_lam2 + 1
  }
  
  if((ci_lam3_l[i] <= true_param[4]) && (ci_lam3_u[i] >= true_param[4])){
    cp_lam3 = cp_lam3 + 1
  }
}

metrics <- c(AL_aph, AL_lam1, AL_lam2, AL_lam3, cp_aph, cp_lam1, cp_lam2, cp_lam3)
metrics <- metrics / rep(n_iter, 8)
print(metrics)
