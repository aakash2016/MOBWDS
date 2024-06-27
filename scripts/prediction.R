# ----------------------------------------------------
# Licensed under the Apache License 2.0.
# Written by Aakash Agrawal (akash.gagrawal@gmail.com)
# ----------------------------------------------------
# %
# Prediction of future failures#

rm(list=ls())
options(warn=-1)

set.seed(43)

source("scripts/utils.R")
library(dplyr)
library(pracma)
library(VGAM)

lifetime_df = read.csv('datasets/device_failure.csv', header = TRUE)
lifetime_df$lifetime <- lifetime_df$lifetime / 150
n = nrow(lifetime_df)

set_num_fails <- function(data) {
  m = table(data$cause)
  
  m_0 <<- 0 # number of failures due to both causes
  m_1 <<- unname(m[1]) # number of failures due to first cause
  m_2 <<- unname(m[2]) # number of failures due to second cause
  m_c <<- unname(m[3]) # number of censored subjects
}

set_num_fails(lifetime_df)

# experimental setup
tau = 3.15
delta = 1.5

# Posterior draws
post_df = read.csv('datasets/posterior_samples_device_failure.csv', header = TRUE)
post_df = post_df[,-1]

SweMinXY <- function(t, params) {
  exp(-params[1]*(t**params[2] + t**params[3] + t**params[4]))
}

cause_ign_rho <- function(params) {
  1 - SweMinXY(tau + delta, params) / SweMinXY(tau, params)
}

jointxy_cause1 <- function(x, y) {
  (x ** (params[2] - 1) * (params[2] * params[1] ** 2) *
     (params[3] * (y ** (params[3] - 1)) + params[4] * (y ** (params[4] - 1))) * 
     exp(-params[1] * (x ** params[2] + y ** params[3] + y ** params[4]))
   )
}

jointxy_cause2 <- function(y, x) {
  (y ** (params[3] - 1) * (params[3] * params[1] ** 2) *
     (params[2] * (x ** (params[2] - 1)) + params[4] * (x ** (params[4] - 1))) * 
     exp(-params[1] * (x ** params[2] + y ** params[3] + x ** params[4]))
  )
}

fxy_1 = function(x) {
  integrate(function(y) {
    (x ** (params[2] - 1) * (params[2] * params[1] ** 2) *
       (params[3] * (y ** (params[3] - 1)) + params[4] * (y ** (params[4] - 1))) * 
       exp(-params[1] * (x ** params[2] + y ** params[3] + y ** params[4]))
    )
  }, x, Inf)$value
}

fxy_2 = function(y) {
  integrate(function(x) {
    (y ** (params[3] - 1) * (params[3] * params[1] ** 2) *
       (params[2] * (x ** (params[2] - 1)) + params[4] * (x ** (params[4] - 1))) * 
       exp(-params[1] * (x ** params[2] + y ** params[3] + x ** params[4]))
    )
  }, y, Inf)$value
}

cause_x_specific_rho <- function(p) {
  params <<- p
  xmin <- tau; xmax <- tau + delta
  ymin <- function(x) x; ymax <- 1000
  # I = integrate(Vectorize(fxy_1), xmin, xmax)$value
  I <- integral2(jointxy_cause1, xmin, xmax, ymin, ymax)$Q
  return(I / SweMinXY(tau, params))
}

cause_y_specific_rho <- function(p) {
  params <<- p
  xmin <- tau; xmax <- tau + delta
  ymin <- function(x) x; ymax <- 1000
  # I = integrate(Vectorize(fxy_2), xmin, xmax)$value
  I <- integral2(jointxy_cause2, xmin, xmax, ymin, ymax)$Q
  return(I / SweMinXY(tau, params))
}

# note: num failures in R -> R + delta is binomial distribution
cause_type = "ignored"
JY_theta <- function(y, p) {
  params <<- p
  if (cause_type == "ignored") {
    rh = cause_ign_rho(params)
  } else if (cause_type == 'cause1') {
    rh = cause_x_specific_rho(params)
  } else {
    rh = cause_y_specific_rho(params)
  }
  pbinom(y, m_c, rh)
}

JY <- function(y) {
  B = 15000
  sum_ = 0
  for (i in 1:B) {
    est = as.numeric(post_df[i,])
    sum_ = sum_ + JY_theta(y, est)
  }
  return(sum_ / B)
}

# log likelihood function
LL <- function(p){
  
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

# MLE
x0 = c(1.40,2.70,1.80,0.74) # initial guess
MLE_optim = constrOptim(x0,
                        LL,
                        grad = NULL,
                        control=list(fnscale=-1),
                        ui=diag(4),
                        ci=rep(0,4)
)

param_mle = MLE_optim$par

## Prediction Frequentist ##
# point estimate
rh_0 = cause_ign_rho(param_mle) # cause ignoring rho
rh_x = cause_x_specific_rho(param_mle) # cause 1 specific rho
rh_y = cause_y_specific_rho(param_mle) # cause 2 specific rho
cat("prediction: ", m_c * c(rh_0, rh_x, rh_y))

# prediction interval
cat("prediction interval (ignoring cause): ", 
    qbinom(0.05, m_c, rh_0), qbinom(0.95, m_c, rh_0))
cat("prediction interval (cause 1): ", 
    qbinom(0.05, m_c, rh_x), qbinom(0.95, m_c, rh_x))
cat("prediction interval (cause 2): ", 
    qbinom(0.05, m_c, rh_y), qbinom(0.95, m_c, rh_y))

## Prediction Bayesian ##
get_bpi <- function(cause_type) {
  cat(cause_type)
  cat('\n')
  
  cause_type <<- cause_type
  bys = JY(0:m_c)
  
  # point estimate
  n_fail = 0
  for (i in 0:m_c) {
    if (i != 0) {
      n_fail = n_fail + (bys[i+1] - bys[i]) * i
    }
  }
  cat("predicted number of failures: ", n_fail)
  cat('\n')
  
  # prediction interval
  ll = max(which((bys < 0.05) == TRUE)) - 1
  if (is.infinite(ll)) {ll = 0}
  ul = min(which((bys > 0.95) == TRUE)) - 1
  cat("prediction interval", ll, ul)
}

get_bpi('ignored')
get_bpi('cause1')
get_bpi('cause2')
