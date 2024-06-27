# ----------------------------------------------------
# Licensed under the Apache License 2.0.
# Written by Aakash Agrawal (akash.gagrawal@gmail.com)
# ----------------------------------------------------
# %
# Hypothesis testing for the data analysis
# LRT -- Likelihood ratio test

rm(list=ls())
options(warn=-1)

source("scripts/utils.R")
library(pracma)
library(nadiv)

lifetime_df = read.csv('datasets/device_failure.csv', header = TRUE)
lifetime_df$lifetime <- lifetime_df$lifetime / 200
n = nrow(lifetime_df)

set_num_fails <- function(data) {
  m = table(data$cause)
  
  m_0 <<- 0 # number of failures due to both causes
  m_1 <<- unname(m[1]) # number of failures due to first cause
  m_2 <<- unname(m[2]) # number of failures due to second cause
  m_c <<- unname(m[3]) # number of censored subjects
}

set_num_fails(lifetime_df)

# Log Likelihood function
# MLE Estimates - LL(H0 U HA) 
LL <- function(x){
  
  sum_ = 0
  sum_j = 0
  
  alpha <<- x[1]
  lam_1 <<- x[2]
  lam_2 <<- x[3]
  lam_3 <<- x[4]
  
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
    (m_0 + m_1 + m_2) * log(alpha) # removed p* 
  
  return(sum_)
}

x0 = c(1.0,1.20,0.88,0.74) # initial guess
MLE_optim = constrOptim(x0,
                        LL,
                        grad = NULL,
                        control=list(fnscale=-1),
                        ui=diag(4),
                        ci=rep(0,4)
)
print(MLE_optim)

# MLE Estimates - LL(H0) 
LL_h0 <- function(x){
  
  sum_ = 0
  sum_j = 0
  
  alpha <<- x[1]
  lam <<- x[2]
  
  for(i in 1:nrow(lifetime_df)) {
    row <- lifetime_df[i, ]
    t = row$lifetime
    
    lambdas = c(lam, lam, lam)
    for (j in lambdas) {
      sum_j = sum_j + log((1 - G(t, j)))
    }
    
    # failure due to the first cause or treatment
    if (row$cause == 1) {
      sum_ = sum_ + log(-R_dash(t, lam, 0) / R(t, lam, 0))
    }
    
    # failure due to the second cause or lack of treatment
    if (row$cause == 2) {
      sum_ = sum_ + log(-R_dash(t, lam, 0) / R(t, lam, 0))
    }
    
    # failure due to both causes
    if (row$cause == 0) {
      sum_ = sum_ + 
        log(-(R_dash(t, lam, 0) / R(t, lam, 0)) - 
              (R_dash(t, lam, lam) / R(t, lam, lam)))
    }
  }
  
  sum_ = sum_ + alpha * sum_j + 
    (m_0 + m_1 + m_2) * log(alpha) # removed p* 
  
  return(sum_)
}

x0 = c(1.0,0.88) # initial guess
MLE_optim_h0 = constrOptim(x0,
                        LL_h0,
                        grad = NULL,
                        control=list(fnscale=-1),
                        ui=diag(2),
                        ci=rep(0,2)
)
print(MLE_optim_h0)

# Under significance level 0.05, 
# we would reject the null hypothesis and conclude that we should use the more complex model.
# compute R factor or the test statistic
tst = -2 * (MLE_optim_h0$value - MLE_optim$value)
pval <- pchisq(tst, df = 1, lower.tail = FALSE)
cat(pval)

# LRT: https://search.r-project.org/CRAN/refmans/nadiv/html/LRTest.html
# https://api.rpubs.com/tomanderson_34/lrt
withBC <- LRTest(full = MLE_optim$value, reduced = MLE_optim_h0$value, df = 1, boundaryCorrection = FALSE)
print(withBC)
