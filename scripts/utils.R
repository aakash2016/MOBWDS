# ----------------------------------------------------
# Licensed under the Apache License 2.0.
# Written by Aakash Agrawal (akash.gagrawal@gmail.com)
# ----------------------------------------------------

# Common Functions
G <- function(t, lambda){
  return(1 - exp(-(t**lambda)))
}

# Distribution function
F <- function(t, alpha, lambda) {
  return(1-(1-G(t, lambda)) ** alpha)
}

#
R <- function(t, lambda_1, lambda_2){
  if (lambda_2 == 0) {
    return(1 - G(t, lambda_1))
  }
  return((1 - G(t, lambda_1)) * (1 - G(t, lambda_2)))
}

# dervative of function G
G_dash <- function(t, lambda){
  return(lambda * (t ** (lambda - 1)) * exp(-(t ** lambda)))
}

# derivative of function R
R_dash <- function(t, lambda_1, lambda_2) {
  # R' = g1 * g2' + g2 * g1' - g1' - g2'
  g1 = G(t, lambda_1)
  g1_d = G_dash(t, lambda_1)
  
  if (lambda_2 == 0) {
    return(-g1_d)
  }
  
  g2 = G(t, lambda_2)
  g2_d = G_dash(t, lambda_2)
  return(g1 * g2_d + g2 * g1_d - g1_d - g2_d)
}

# Num fails for each cause
set_num_fails <- function(data) {
  m = table(data$cause)
  
  m_0 <<- unname(m[1]) # number of failures due to both causes
  m_1 <<- unname(m[2]) # number of failures due to first cause
  m_2 <<- unname(m[3]) # number of failures due to second cause
  m_c <<- unname(m[4]) # number of censored subjects
}

# CI from CV, params
get_CI_OFIM <- function(params, CV) {
  ci_aph_l = params[1] - 1.96 * sqrt(CV[1, 1])
  ci_aph_u = params[1] + 1.96 * sqrt(CV[1, 1])
  ci_lam1_l = params[2] - 1.96 * sqrt(CV[2, 2])
  ci_lam1_u = params[2] + 1.96 * sqrt(CV[2, 2])
  ci_lam2_l = params[3] - 1.96 * sqrt(CV[3, 3])
  ci_lam2_u = params[3] + 1.96 * sqrt(CV[3, 3])
  ci_lam3_l = params[4] - 1.96 * sqrt(CV[4, 4])
  ci_lam3_u = params[4] + 1.96 * sqrt(CV[4, 4])  
  return(c(ci_aph_l, ci_aph_u, ci_lam1_l, ci_lam1_u, 
           ci_lam2_l, ci_lam2_u, ci_lam3_l, ci_lam3_u))
}

# Save the results in csv file
write_data <- function(n_iter, 
                       MLE_est, 
                       alpha_est, lam_1_est, lam_2_est, lam_3_est, 
                       ci_aph_l, ci_aph_u, ci_lam1_l, ci_lam1_u,
                       ci_lam2_l, ci_lam2_u, ci_lam3_l, ci_lam3_u,
                       convergence) {
  cols = c("MLE", 
           "alpha", "lambda_1", "lambda_2", "lambda_3", 
           "CI_alpha_l", "CI_alpha_u", "CI_lambda1_l", "CI_lambda1_u",
           "CI_lambda2_l", "CI_lambda2_u", "CI_lambda3_l", "CI_lambda3_u",
           "convergence")
  
  df = data.frame(matrix(nrow = n_iter, ncol = length(cols)))
  colnames(df) = cols
  
  df$MLE = MLE_est
  df$alpha = alpha_est
  df$lambda_1 = lam_1_est
  df$lambda_2 = lam_2_est
  df$lambda_3 = lam_3_est
  df$CI_alpha_l = ci_aph_l
  df$CI_alpha_u = ci_aph_u
  df$CI_lambda1_l = ci_lam1_l
  df$CI_lambda1_u = ci_lam1_u
  df$CI_lambda2_l = ci_lam2_l
  df$CI_lambda2_u = ci_lam2_u
  df$CI_lambda3_l = ci_lam3_l
  df$CI_lambda3_u = ci_lam3_u
  df$convergence = convergence
  
  return(df)
}
