# ----------------------------------------------------
# Licensed under the Apache License 2.0.
# Written by Aakash Agrawal (akash.gagrawal@gmail.com)
# ----------------------------------------------------
# %
# Simulation Study for Sensitivity Analysis #

rm(list=ls())
options(warn=-1)

source("scripts/utils.R")
library(VGAM)

n_dpts = 200 # sample size
censoring = FALSE # whether to have censoring or not
tau = 0.5 # censoring factor

n_iter = 1000 # Number of monte carlo replications
perc_25 = 0.1 # found using 10k samples generated
perc_75 = 0.33 
t_range = seq(perc_25, perc_75, by=0.25*(perc_75-perc_25))
true_param = c(2.35,1.11,1.92,1.63)

# unit reliability
relative_ur_M1 = relative_ur_M2 = matrix(nrow=n_iter,ncol=5)

# Estimates
alpha_est = lam_1_est = lam_2_est = lam_3_est = MLE_est = rep(0,n_iter)

# Convergence of MLEs
convergence = rep(0,n_iter)

# relaibilty of a unit
unit_reliability <- function(t, params) {
  return(exp(-params[1] * ((t ** params[2]) + (t ** params[3]) + (t ** params[4]))))
}
true_ur =  unit_reliability(t_range, true_param)

# functions to get the probability factors
f_integral_weibull_1 <- function(t) {
  return(alpha * lam_1 * (t ** (lam_1 - 1)) * 
           exp(-alpha * ((t ** lam_1) + (t ** lam_2) + (t ** lam_3))))
}

f_integral_weibull_2 <- function(t) {
  return(alpha * lam_2 * (t ** (lam_2 - 1)) * 
           exp(-alpha * ((t ** lam_1) + (t ** lam_2) + (t ** lam_3))))
}

f_integral_gompertz_1 <- function(t) {
  return(alpha * lam_1 * exp((lam_1 * t) + 
                               alpha * (3 - exp(lam_1 * t) - exp(lam_2 * t) - exp(lam_3 * t))))
}

f_integral_gompertz_2 <- function(t) {
  return(alpha * lam_2 * exp((lam_2 * t) + 
                               alpha * (3 - exp(lam_1 * t) - exp(lam_2 * t) - exp(lam_3 * t))))
}

f_integral_lomax_1 <- function(t) {
  return(alpha * lam_1 * ((1 + lam_1 * t) ** -(alpha + 1)) * 
           ((1 + lam_2 * t) ** -alpha) * ((1 + lam_3 * t) ** -alpha))
}

f_integral_lomax_2 <- function(t) {
  return(alpha * lam_2 * ((1 + lam_2 * t) ** -(alpha + 1)) * 
           ((1 + lam_1 * t) ** -alpha) * ((1 + lam_3 * t) ** -alpha))
}

# Log Likelihood function for M1: MOBW
LL_M1 <- function(x){
  
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

# Log Likelihood function M2: Independent We
LL_M2 <- function(x){
  
  sum_1 = 0
  sum_2 = 0
  sum_j = 0
  
  aph_1 <- x[1]
  aph_2 <- x[2]
  lam <- x[3]
  
  for(i in 1:nrow(lifetime_df)) {
    row <- lifetime_df[i, ]
    t = row$lifetime
    
    sum_j = sum_j + t ^ aph_1 + t ^ aph_2
    
    if (row$cause == 1) {
      sum_1 = sum_1 + log(t)
    }
    
    if (row$cause == 2) {
      sum_2 = sum_2 + log(t)
    }
  }
  
  sum_ = (aph_1 - 1) * sum_1 + 
    (aph_2 - 1) * sum_2 - 
    lam * sum_j + 
    m_1 * log(aph_1) +
    m_2 * log(aph_2) +
    (m_1 + m_2) * log(lam)
  
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
  
  # MLE: M1
  MLE_M1 = try(
    {constrOptim(true_param,
                 LL_M1,
                 grad=NULL,
                 control=list(fnscale=-1),
                 ui=diag(4),
                 ci=rep(0,4))},
    silent = TRUE
  )
  
  # modify dataset and randomly split tites 
  lifetime_df[lifetime_df$cause == 0, ]$cause = (
    1 + as.integer(runif(length(lifetime_df[lifetime_df$cause == 0, ]$cause)) > 0.5)
    )
  
  lifetime_df <<- lifetime_df
  m = table(lifetime_df$cause)
  m_0 <<- 0
  m_1 <<- unname(m[1])
  m_2 <<- unname(m[2])
  print(c(m_0, m_1, m_2, m_c))
  
  # MLE: M2
  MLE_M2 = try(
    {constrOptim(c(true_param[2], true_param[3], true_param[1]),
                 LL_M2,
                 grad = NULL,
                 control=list(fnscale=-1),
                 ui=diag(3),
                 ci=rep(0,3))}, 
    silent = TRUE
  )
  
  # if error then skip
  if(inherits(MLE_M1, 'try-error') | inherits(MLE_M2, 'try-error')) {
    mle_error = mle_error + 1
    next
  }
  
  # Divergence Cases
  if (abs(MLE_M1$par[4] - true_param[4]) > 10) {
    next
  }
  
  # log the values
  MLE_est[i] = MLE_M1$value
  alpha_est[i] = aph = MLE_M1$par[1]
  lam_1_est[i] = l1 = MLE_M1$par[2]
  lam_2_est[i] = l2 = MLE_M1$par[3]
  lam_3_est[i] = l3 = MLE_M1$par[4]
  
  # t calculation
  # lms = unname(quantile(lifetime_df$lifetime, c(0.25, 0.75)))
  ur_M1 = unit_reliability(t_range, MLE_M1$par)
  ur_M2 = unit_reliability(t_range, c(MLE_M2$par[3], MLE_M2$par[1], MLE_M2$par[2], 0))
  
  relative_ur_M1[i,] = abs(true_ur - ur_M1) / true_ur
  relative_ur_M2[i,] = abs(true_ur - ur_M2) / true_ur
  
  i = i + 1
}

end.time <- Sys.time()
time.taken <- end.time - start.time
cat("time taken for simulation: ", time.taken)

cat(t_range)
cat(colMeans(relative_ur_M1))
cat(colMeans(relative_ur_M2))
