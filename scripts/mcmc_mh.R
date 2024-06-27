# ----------------------------------------------------
# Licensed under the Apache License 2.0.
# Written by Aakash Agrawal (akash.gagrawal@gmail.com)
# ----------------------------------------------------
# %
# Markov Chain Monte Carlo #

rm(list=ls())
options(warn=-1)

set.seed(43)

source("scripts/utils.R")
library(dplyr)
library(VGAM)

# Electric Surge dataset
lifetime_df = read.csv('datasets/device_failure.csv', header = TRUE)
lifetime_df$lifetime <- lifetime_df$lifetime / 150
set_num_fails(lifetime_df)

# log likelihood function
loglikelihood <- function(lifetime_df, p){
  
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
logprior <- function(p) {
  c1 = c2 = a = b = 0.005
  alpha = rep(1.2, 3)
  log(GA_lambda(p[1], c1, c2) * GD_alpha(p[2:4], alpha, a, b))
}

# posterior density log scale
posterior <- function(x, p) {
  loglikelihood(x, p) + logprior(p)
}

# var cov matrix
rho = 0.1
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
  i = 1
  samples = x_curr
  while (i < N) {
    # track progress
    if(i %% 100 == 0){
      cat("iteration num: ", i)
      cat("\n")
    }
    
    # propose new candidate value using folded normal
    theta_star <- rmvfoldednormal(x_curr, COV)
    
    # calculate the ratio R
    pstar <- try(
      {
        posterior(lifetime_df, theta_star)
      },
      silent = TRUE
    )
    
    pprev <- try(
      {
        posterior(lifetime_df, x_curr) 
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
    
    i = i + 1
  }
  
  return(samples)
}

# Experimental Setup
n_steps = 15000 # markov steps
burn_steps = 500
N = burn_steps + n_steps
x_curr = c(1.40,0.95,1.85,0.60) # initial guess
samples = MCMC(x_curr, N) # mcmc-mh sampling
df = data.frame(samples) %>% slice(burn_steps+1:N)
row.names(df) <- NULL

# save posterior samples
## write.csv(x=df, file = 'datasets/posterior_samples_device_failure.csv')

# bayes estimates
lam_est <- mean(df[,1])
aph1_est <- mean(df[,2])
aph2_est <- mean(df[,3])
aph0_est <- mean(df[,4])

# posterior variance
var_lam = mean((lam_est - df[,1])**2)
var_aph1 = mean((aph1_est - df[,2])**2)
var_aph2 = mean((aph2_est - df[,3])**2)
var_aph0 = mean((aph0_est - df[,4])**2)

# 95% credible intervals
credible_ci = function(x){
  ci = unname(quantile(x, c(0.05, 0.95), na.rm = FALSE))
  return(ci)
}

lam_ci <- credible_ci(df[,1])
aph1_ci <- credible_ci(df[,2])
aph2_ci <- credible_ci(df[,3])
aph0_ci <- credible_ci(df[,4])

# trace plots
par(mfrow=c(2,2))
plot(1:n_steps, df[,4], type="l", xlab = "cycles", ylab = expression(alpha[0]))
abline(h=aph0_est, col="blue", lty=2)

plot(1:n_steps, df[,2], type="l", xlab = "cycles", ylab = expression(alpha[1]))
abline(h=aph1_est, col="blue", lty=2)

plot(1:n_steps, df[,3], type="l", xlab = "cycles", ylab = expression(alpha[2]))
abline(h=aph2_est, col="blue", lty=2)

plot(1:n_steps, df[,1], type="l", xlab = "cycles", ylab = expression(lambda))
abline(h=lam_est, col="blue", lty=2)

# Autocorrelation Plots
acf(df[,4], type = 'correlation', plot=TRUE, lag.max = 200, main=expression(alpha[0]))
acf(df[,2], type = 'correlation', plot=TRUE, lag.max = 200, main=expression(alpha[1]))
acf(df[,3], type = 'correlation', plot=TRUE, lag.max = 200, main=expression(alpha[2]))
acf(df[,1], type = 'correlation', plot=TRUE, lag.max = 200, main=expression(lambda))

# Gelman Measure
x1 = c(1.92,0.81,1.09,1.33)
df_1 = data.frame(MCMC(x1, N)) %>% slice(burn_steps+1:N)
row.names(df_1) <- NULL

x2 = c(0.88,0.66,1.24,0.75)
df_2 = data.frame(MCMC(x2, N)) %>% slice(burn_steps+1:N)
row.names(df_2) <- NULL

B_lam = sd(c(df_1[,1], df_2[,1])) # between chain variability
W_lam = 0.5 * (sd(df_1[,1]) + sd(df_2[,1])) # within chain variability
R_lam = B_lam / W_lam

B_aph1 = sd(c(df_1[,2], df_2[,2]))
W_aph1 = 0.5 * (sd(df_1[,2]) + sd(df_2[,2])) 
R_aph1 = B_aph1 / W_aph1

B_aph2 = sd(c(df_1[,3], df_2[,3]))
W_aph2 = 0.5 * (sd(df_1[,3]) + sd(df_2[,3]))
R_aph2 = B_aph2 / W_aph2

B_aph0 = sd(c(df_1[,4], df_2[,4]))
W_aph0 = 0.5 * (sd(df_1[,4]) + sd(df_2[,4]))
R_aph0 = B_aph0 / W_aph0

cat(R_lam, R_aph1, R_aph2, R_aph0)

