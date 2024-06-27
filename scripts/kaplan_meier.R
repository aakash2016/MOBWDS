# ----------------------------------------------------
# Licensed under the Apache License 2.0.
# Written by Aakash Agrawal (akash.gagrawal@gmail.com)
# ----------------------------------------------------
# %
# data analysis - MTTF and Kaplan Meier #

rm(list=ls())
options(warn=-1)

library(dplyr)
library(survival)
library(VGAM)

source("scripts/utils.R")

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
lifetime_df$ind = 1  
lifetime_df$ind[lifetime_df$cause == 3] = 0

# check for bimodality in the histogram
hist(lifetime_df$lifetime, nclass=50)

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

# MLE Estimates
x0 = c(1.40,2.70,1.80,0.74) # initial guess
MLE_optim = constrOptim(x0,
                        LL,
                        grad = NULL,
                        control=list(fnscale=-1),
                        ui=diag(4),
                        ci=rep(0,4)
)

params = MLE_optim$par

# Bayes Estimates
post_df = read.csv('datasets/posterior_samples_device_failure.csv', header = TRUE)
post_df = post_df[,-1]

lam_est <- mean(post_df[,1])
aph1_est <- mean(post_df[,2])
aph2_est <- mean(post_df[,3])
aph0_est <- mean(post_df[,4])

bayes_par = c(lam_est, aph1_est, aph2_est, aph0_est)

## Parametric Estimate of Survival (weibull)
SweMinXY <- function(t, params) {
  return(exp(-params[1]*(t**params[2] + t**params[3] + t**params[4])))
}

# Kaplan-Meier fit -- Non Parametric
KMfit_MinXY = survfit(Surv(lifetime, ind) ~ 1, se.fit=FALSE, data=lifetime_df)

par(mfrow=c(1,1))
t = seq(0.0, 3.15, 0.001)
plot(KMfit_MinXY,
     mgp=c(2.25, 1, 0.),
     xlab = expression("T (units normalized)"),
     ylab = "Estimated Survival Probability")
lines(t, SweMinXY(t, params), type = 'l', col="blue")
lines(t, SweMinXY(t, bayes_par), type = 'l', col="red")
legend(x = "bottomleft", lty = 1,
       col= c("blue","red", "black"),
       legend=c("Frequentist", "Bayesian", "Kaplan-Meier"),
       box.lty=0,
       box.lwd=0,
       cex=0.8)

# dev.copy(jpeg, filename="plots/survest_minxy.jpg")
# dev.off()

# MTTF Frquentist
mttf_freq = integrate(SweMinXY, lower = 0, upper = Inf, params = params)$value * 150

# MTTF Bayesian
N = 15000
mttf_mcmc = rep(0, N)
for (i in 1:N) {
  est = as.numeric(post_df[i,])
  mttf_mcmc[i] = integrate(SweMinXY, lower = 0, upper = Inf, params = est)$value
}

mttf_bayes = mean(mttf_mcmc) * 150
cat(mttf_freq, mttf_bayes)
