# ----------------------------------------------------
# Licensed under the Apache License 2.0.
# Written by Aakash Agrawal (akash.gagrawal@gmail.com)
# ----------------------------------------------------
# %
# Bayesian Simulation #

rm(list=ls())
options(warn=-1)

source("scripts/utils.R")
library(dplyr)
library(pracma)
library(VGAM)

n_dpts = 400 # sample size
censoring = TRUE # whether to have censoring or not
tau = 0.3 # censoring factor

n_steps = 1000 # markov steps
burn_steps = 50
N = burn_steps + n_steps

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

# log likelihood function
loglikelihood <- function(p){
  
  ## alpha, lambda_1, lambda_2, lambda_3 are the estimates we want to find
  sum_ = 0
  sum_j = 0
  
  alpha <- p[1]
  lam_1 <- p[2]
  lam_2 <- p[3]
  lam_3 <- p[4]
  
  # functions to get the probability factors
  f_integral_1 <- function(t) {
    return(alpha * lam_1 * (t ** (lam_1 - 1)) * 
             exp(-alpha * ((t ** lam_1) + (t ** lam_2) + (t ** lam_3))))
  }
  
  f_integral_2 <- function(t) {
    return(alpha * lam_2 * (t ** (lam_2 - 1)) * 
             exp(-alpha * ((t ** lam_1) + (t ** lam_2) + (t ** lam_3))))
  }
  
  p1 = integrate(f_integral_1, lower = 0, upper = Inf)$value
  p2 = integrate(f_integral_2, lower = 0, upper = Inf)$value
  
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

# gamma-dirichlet prior for different shape's alpha
GD_alpha <- function(theta, alpha, a, b) {
  theta_n <- theta / sum(theta)
  A <- sum(alpha) 
  return (
    (theta_n[1] ** (alpha[1] - 1))   # aph1
    * (theta_n[2] ** (alpha[2] - 1)) # aph2
    * (theta_n[3] ** (alpha[3] - 1)) # aph0
    * (exp(-b * A)) 
    * ((A) ** (a - 3))
  )
}

# gamma prior for common scale lambda
GA_lambda <- function(lambda, c1, c2) {
  exp(-c1 * lambda) * (lambda ** (c2 - 1))
}

# prior density
# informative GD prior for joint alpha's
kt = 2.0
aps = kt * true_param[2:4]
ap_0 = sum(aps)
var_aps = (aps * (ap_0 - aps)) / ((ap_0 + 1) * ap_0 ** 2) # variance
cat("variance of GD hyperparams: ", var_aps)

logprior <- function(p) {
  c1 = c2 = a = b = 0.005 # GA lambda
  log(GA_lambda(p[1], c1, c2) * GD_alpha(p[2:4], aps, a, b))
}

# posterior density log scale
posterior <- function(p) {
  loglikelihood(p) + logprior(p)
}

# var cov matrix
rho = 0.5
COV = matrix( 
  c(
    rho, 0., 0., 0.,
    0., rho, 0., 0.,
    0., 0., rho, 0.,
    0., 0., 0., rho
  ),
  nrow = 4,
  ncol = 4,
  byrow = TRUE
)

# density of abs normal
dabsnorm <- function(y, mu, sig) {
  if (y >= 0) {
    (exp(-0.5 * (y - mu) ** 2 / sig ** 2) + 
       exp(-0.5 * (-y - mu) ** 2 / sig ** 2)) / (sqrt(2 * pi) * sig) 
  } else {
    0
  }
}

# multivariate density of abs norm
mvdabsnorm <- function(x, MU, S) {
  sum = 0
  for (i in 1:4) {
    sum = sum + log(dabsnorm(x[i], MU[i], S[i,i]))
  }
  return(exp(sum))
}

# sample from multivariate abs normal
rmvfoldednormal <- function(x_curr, COV) {
  p1 <- rfoldnorm(1, x_curr[1], COV[1,1])
  p2 <- rfoldnorm(1, x_curr[2], COV[2,2])
  p3 <- rfoldnorm(1, x_curr[3], COV[3,3])
  p4 <- rfoldnorm(1, x_curr[4], COV[4,4])
  return(c(p1, p2, p3, p4))
}

# MCMC - MH
MCMC <- function(x_curr, N) {
  j = 1
  samples = x_curr
  while (j < N) {
    # propose new candidate value using folded normal
    theta_star <- rmvfoldednormal(x_curr, COV)
    
    # calculate the ratio R
    pstar <- try(
      {
        posterior(theta_star)
      },
      silent = TRUE
    )
    
    pprev <- try(
      {
        posterior(x_curr) 
      }, 
      silent = TRUE
    ) 
    
    # if error then skip
    if(inherits(pstar, 'try-error') || 
       inherits(pprev, 'try-error') || 
       is.nan(pstar)) {
      next
    }
    
    logR <- pstar - pprev + (
      log(mvdabsnorm(x_curr, theta_star, COV)) - 
        log(mvdabsnorm(theta_star, x_curr, COV))
    )
    
    RA <- min(1, exp(logR))
    
    # accept or reject
    if (runif(1) < RA) {
      x_curr <- theta_star
    }
    
    samples <- rbind(samples, x_curr)
    
    j = j + 1
  }
  
  return(samples)
}

i = 1
hessian_error = 0
start.time <- Sys.time()
while (i < n_iter) {
  # track progress
  if(i %% 1 == 0){
    cat("iteration num: ", i)
    cat("\n")
  }
  
  # dataset
  df = get_dataset(true_param)
  print(c(m_0, m_1, m_2, m_c))
  
  lifetime_df <<- df # make df global
  
  # MAE
  samples = MCMC(true_param, N) # mcmc-mh sampling
  MC_df = data.frame(samples) %>% slice(burn_steps+1:N)
  row.names(MC_df) <- NULL
  
  # bayes estimates
  alpha_est[i] = aph = mean(MC_df[,1])
  lam_1_est[i] = l1 = mean(MC_df[,2])
  lam_2_est[i] = l2 = mean(MC_df[,3])
  lam_3_est[i] = l3 = mean(MC_df[,4])
  
  cat("alpha_0", l3 - true_param[4])
  cat("\n")
  
  # covariance matrix is inverse of negative of hessian
  CV = inv(-hessian(loglikelihood, c(aph, l1, l2, l3)))
  
  # if diagonal element less than 0
  if (any(diag(CV) < 0) | all(is.na(CV))) {
    hessian_error = hessian_error + 1
    next
  }
  
  # Divergence Cases
  if (abs( l3 - true_param[4]) > 10) {
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
