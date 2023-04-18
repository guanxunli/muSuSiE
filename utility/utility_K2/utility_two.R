#### Define functions
## get sigma0
lBF_model_two <- function(lsigma02, prior_pi, z2_1, s2_1, z2_2, s2_2) {
  sigma02 <- exp(lsigma02)
  # data set 1
  tmp1_1 <- log(sqrt(s2_1 / (sigma02 + s2_1)))
  tmp2_1 <- z2_1 / 2 * sigma02 / (sigma02 + s2_1)
  lBF_1 <- tmp1_1 + tmp2_1
  # data set 2
  tmp1_2 <- log(sqrt(s2_2 / (sigma02 + s2_2)))
  tmp2_2 <- z2_2 / 2 * sigma02 / (sigma02 + s2_2)
  lBF_2 <- tmp1_2 + tmp2_2
  # combine
  lBF <- c(lBF_1, lBF_2, lBF_1 + lBF_2, 0)
  maxlBF <- max(lBF)
  wBF <- exp(lBF - maxlBF)
  wBF_sum <- sum(prior_pi * wBF)
  return(-maxlBF - log(wBF_sum))
}

sigma0_opt_two <- function(lsigma02_int, prior_pi, z2_1, s2_1, z2_2, s2_2, b_hat_1, b_hat_2) {
  tmp1 <- lBF_model_two(
    lsigma02 = lsigma02_int, prior_pi = prior_pi, z2_1 = z2_1,
    s2_1 = s2_1, z2_2 = z2_2, s2_2 = s2_2
  )
  # lsigma02 <- optim(
  #   par = log(max(c(b_hat_1^2 - s2_1, 1, b_hat_2^2 - s2_2))), fn = lBF_model_two,
  #   method = "Brent", lower = -30, upper = 15, prior_pi = prior_pi, z2_1 = z2_1,
  #   s2_1 = s2_1, z2_2 = z2_2, s2_2 = s2_2
  # )$par
  lsigma02 <- optimize(
    lBF_model_two,
    lower = -30, upper = 15, prior_pi = prior_pi, z2_1 = z2_1,
    s2_1 = s2_1, z2_2 = z2_2, s2_2 = s2_2
  )$minimum
  tmp2 <- lBF_model_two(
    lsigma02 = lsigma02, prior_pi = prior_pi, z2_1 = z2_1,
    s2_1 = s2_1, z2_2 = z2_2, s2_2 = s2_2
  )
  if (tmp2 < tmp1) {
    return(exp(lsigma02))
  } else {
    return(exp(lsigma02_int))
  }
}

## Calculate the -KL divergence
KL_fun_two <- function(X_scale_1, X_scale2_1, Y_1, X_scale_2, Y_2, X_scale2_2,
                       sigma2, b_1, b2_1, b_2, b2_2, lBF) {
  n1 <- length(Y_1)
  n2 <- length(Y_2)
  tmp1_1 <- sum(dnorm(Y_1, mean = 0, sd = sqrt(sigma2), log = TRUE))
  tmp1_2 <- sum(dnorm(Y_2, mean = 0, sd = sqrt(sigma2), log = TRUE))
  tmp3_1 <- n1 / 2 * log(2 * pi * sigma2)
  tmp3_2 <- n2 / 2 * log(2 * pi * sigma2)
  tmp4_1 <- 1 / (2 * sigma2) * (crossprod(Y_1) - 2 * crossprod(Y_1, X_scale_1 %*% b_1) + sum(X_scale2_1 %*% b2_1))
  tmp4_2 <- 1 / (2 * sigma2) * (crossprod(Y_2) - 2 * crossprod(Y_2, X_scale_2 %*% b_2) + sum(X_scale2_2 %*% b2_2))
  return(tmp1_1 + tmp1_2 + lBF + tmp3_1 + tmp3_2 + tmp4_1 + tmp4_2)
}

KL_fun_two_graph <- function(XtY_1, XtXbeta_use_1, X_scale2_1,
                             XtY_2, XtXbeta_use_2, X_scale2_2,
                             sigma2, b_1, b2_1, b_2, b2_2, lBF) {
  tmp1_1 <- -2 * (crossprod(b_1, XtY_1) - crossprod(b_1, XtXbeta_use_1))
  tmp1_2 <- -2 * (crossprod(b_2, XtY_2) - crossprod(b_2, XtXbeta_use_2))
  tmp2_1 <- sum(X_scale2_1 %*% b2_1)
  tmp2_2 <- sum(X_scale2_2 %*% b2_2)
  return(lBF + 1 / (2 * sigma2) * (tmp1_1 + tmp1_2 + tmp2_1 + tmp2_2))
}

## Calculate ERSS
ERSS_fun_single <- function(X_scale, X_scale2, Y, b_mat, b2_mat) {
  mu_lmat <- X_scale %*% b_mat
  mu2_lmat <- X_scale2 %*% b2_mat
  mu_pred <- rowSums(mu_lmat)
  res_tmp <- sum((Y - mu_pred)^2)
  var_sum <- sum(mu2_lmat - mu_lmat^2)
  return(res_tmp + var_sum)
}