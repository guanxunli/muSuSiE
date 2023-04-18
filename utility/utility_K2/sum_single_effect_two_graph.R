source("utility/utility_K2/utility_two.R")
sum_single_effect_two <- function(X_1, Y_1, X_2, Y_2, scale_x = TRUE, intercept = TRUE,
                                  sigma02_int = NULL, sigma2_int = NULL, prior_vec = NULL,
                                  L = NULL, itermax = 100, tol = 1e-4, sigma0_low_bd = 1e-8,
                                  residual_variance_lowerbound = NULL) {
  ## Initialization
  p <- ncol(X_1)
  if (p != ncol(X_2)) stop("The number of features should be same!")
  n1 <- nrow(X_1)
  n2 <- nrow(X_2)

  # Initialize sigma
  if (is.null(sigma2_int)) sigma2_int <- as.numeric(var(c(Y_1, Y_2)))
  if (is.null(sigma02_int)) sigma02_int <- 0.2 * sigma2_int
  if (is.null(L)) L <- min(10, p)
  if (is.null(residual_variance_lowerbound)) residual_variance_lowerbound <- 1e-4

  ## data preprocess
  X_scale_1 <- scale(X_1, center = intercept, scale = scale_x)
  X_scale_2 <- scale(X_2, center = intercept, scale = scale_x)

  if (intercept) {
    mean_Y_1 <- mean(Y_1)
    mean_Y_2 <- mean(Y_2)
  } else {
    mean_Y_1 <- 0
    mean_Y_2 <- 0
  }
  Y_1 <- Y_1 - mean_Y_1
  Y_2 <- Y_2 - mean_Y_2

  ## data set 1
  XtX_1 <- crossprod(X_scale_1)
  X_scale2_1 <- X_scale_1 * X_scale_1
  X2_1 <- colSums(X_scale2_1)
  XtY_1 <- crossprod(X_scale_1, Y_1)
  # data set 2
  XtX_2 <- crossprod(X_scale_2)
  X_scale2_2 <- X_scale_2 * X_scale_2
  X2_2 <- colSums(X_scale2_2)
  XtY_2 <- crossprod(X_scale_2, Y_2)

  # Initialize prior
  if (is.null(prior_vec)) {
    prior_vec <- c(1 / (2 * p^1.5), 1 / (p^2))
  }
  prior_pi <- c(rep(prior_vec[1], 2 * p), rep(prior_vec[2], p))
  prior_pi <- c(prior_pi, 1 - sum(prior_pi))

  # initialize ELBO
  ELBO <- rep(NA, itermax + 1)
  ELBO[1] <- -Inf
  sigma2 <- sigma2_int
  sigma02_vec <- rep(sigma02_int, L)

  # Save matrix
  b_mat_1 <- matrix(0, nrow = p, ncol = L)
  b2_mat_1 <- matrix(0, nrow = p, ncol = L)
  b_mat_2 <- matrix(0, nrow = p, ncol = L)
  b2_mat_2 <- matrix(0, nrow = p, ncol = L)
  alpha_mat_1 <- matrix(0, nrow = p, ncol = L)
  alpha_mat_2 <- matrix(0, nrow = p, ncol = L)
  alpha_mat <- matrix(0, nrow = 3 * p, ncol = L)

  # Begin iteration
  beta_hat_1 <- rowSums(b_mat_1)
  beta_hat_2 <- rowSums(b_mat_2)
  for (iter in seq_len(itermax)) {
    KL_div <- 0
    for (l in seq_len(L)) {
      ## data set 1
      beta_use_1 <- beta_hat_1 - b_mat_1[, l]
      # update parameters
      XtXbeta_use_1 <- XtX_1 %*% beta_use_1
      XtYtmp_1 <- XtY_1 - XtXbeta_use_1
      b_hat_1 <- XtYtmp_1 / X2_1
      s2_1 <- sigma2 / X2_1
      z2_1 <- b_hat_1^2 / s2_1
      ## data set 2
      beta_use_2 <- beta_hat_2 - b_mat_2[, l]
      # update parameters
      XtXbeta_use_2 <- XtX_2 %*% beta_use_2
      XtYtmp_2 <- XtY_2 - XtXbeta_use_2
      b_hat_2 <- XtYtmp_2 / X2_2
      s2_2 <- sigma2 / X2_2
      z2_2 <- b_hat_2^2 / s2_2
      # calculate sigma0
      lsigma02_int <- max(log(sigma02_vec[l]), -30)
      sigma02 <- sigma0_opt_two(lsigma02_int, prior_pi, z2_1, s2_1, z2_2, s2_2, b_hat_1, b_hat_2)
      sigma02_vec[l] <- sigma02
      ## Get Bayesian Factor
      # data set 1
      tmp1_1 <- log(sqrt(s2_1 / (sigma02 + s2_1)))
      tmp2_1 <- z2_1 / 2 * sigma02 / (sigma02 + s2_1)
      lBF_1 <- tmp1_1 + tmp2_1
      # data set 2
      tmp1_2 <- log(sqrt(s2_2 / (sigma02 + s2_2)))
      tmp2_2 <- z2_2 / 2 * sigma02 / (sigma02 + s2_2)
      lBF_2 <- tmp1_2 + tmp2_2
      # get bayesian factor
      lBF <- c(lBF_1, lBF_2, lBF_1 + lBF_2, 0)
      maxlBF <- max(lBF)
      wBF <- exp(lBF - maxlBF)
      wBF_sum <- sum(prior_pi * wBF)
      lBF_model <- maxlBF + log(wBF_sum)
      ## Get posterior
      post_alpha <- prior_pi * wBF / wBF_sum
      alpha_mat[, l] <- post_alpha[-length(post_alpha)]
      # data set 1
      post_sigma2_1 <- 1 / (1 / s2_1 + 1 / sigma02)
      post_mu_1 <- post_sigma2_1 / s2_1 * b_hat_1
      # data set 2
      post_sigma2_2 <- 1 / (1 / s2_2 + 1 / sigma02)
      post_mu_2 <- post_sigma2_2 / s2_2 * b_hat_2
      ## Calculate posterior mean
      # data set 1
      alpha_mat_1[, l] <- post_alpha[1:p] + post_alpha[(2 * p + 1):(3 * p)]
      b_mat_1[, l] <- alpha_mat_1[, l] * post_mu_1
      b2_mat_1[, l] <- alpha_mat_1[, l] * (post_mu_1^2 + post_sigma2_1)
      # data set 2
      alpha_mat_2[, l] <- post_alpha[(p + 1):(2 * p)] + post_alpha[(2 * p + 1):(3 * p)]
      b_mat_2[, l] <- alpha_mat_2[, l] * post_mu_2
      b2_mat_2[, l] <- alpha_mat_2[, l] * (post_mu_2^2 + post_sigma2_2)
      ## calculate the -KL divergence
      KL_div <- KL_div + KL_fun_two_graph(
        XtY_1 = XtY_1, XtXbeta_use_1 = XtXbeta_use_1, X_scale2_1 = X_scale2_1,
        XtY_2 = XtY_2, XtXbeta_use_2 = XtXbeta_use_2, X_scale2_2 = X_scale2_2,
        sigma2 = sigma2, b_1 = b_mat_1[, l], b2_1 = b2_mat_1[, l],
        b_2 = b_mat_2[, l], b2_2 = b2_mat_2[, l], lBF = lBF_model
      )
      beta_hat_1 <- beta_use_1 + b_mat_1[, l]
      beta_hat_2 <- beta_use_2 + b_mat_2[, l]
    }
    # calculate ELBO
    ERSS_1 <- ERSS_fun_single(X_scale = X_scale_1, X_scale2 = X_scale2_1, Y = Y_1, b_mat = b_mat_1, b2_mat = b2_mat_1)
    ERSS_2 <- ERSS_fun_single(X_scale = X_scale_2, X_scale2 = X_scale2_2, Y = Y_2, b_mat = b_mat_2, b2_mat = b2_mat_2)
    ERSS <- ERSS_1 + ERSS_2
    ELBO[iter + 1] <- -(n1 + n2) / 2 * log(2 * pi * sigma2) - 1 / (2 * sigma2) * ERSS + KL_div
    # estimate sigma2
    sigma2 <- max(ERSS / (n1 + n2), residual_variance_lowerbound)
    if (ELBO[iter + 1] - ELBO[iter] < 1e-4) break
  }
  ELBO <- as.numeric(na.omit(ELBO[-1]))
  # select effect index
  index_L <- which(sigma02_vec > sigma0_low_bd)
  ## return results
  res <- list()
  res$ELBO <- ELBO
  res$sigma2 <- sigma2
  res$sigma02_vec <- sigma02_vec

  if (length(index_L) > 0) {
    # all alpha
    res$alpha <- 1 - matrixStats::rowProds(1 - alpha_mat[, index_L, drop = FALSE])
    # data set 1
    res$alpha_1 <- 1 - matrixStats::rowProds(1 - alpha_mat_1[, index_L, drop = FALSE])
    res$post_mean1 <- rowSums(b_mat_1[, index_L, drop = FALSE])
    res$Xb_1 <- mean_Y_1 + X_scale_1 %*% res$post_mean1
    # data set 2
    res$alpha_2 <- 1 - matrixStats::rowProds(1 - alpha_mat_2[, index_L, drop = FALSE])
    res$post_mean2 <- rowSums(b_mat_2[, index_L, drop = FALSE])
    res$Xb_2 <- mean_Y_2 + X_scale_2 %*% res$post_mean2
  } else {
    # all alpha
    res$alpha <- rep(0, 3 * p)
    # data set 1
    res$alpha_1 <- rep(0, p)
    res$post_mean1 <- rep(0, p)
    res$Xb_1 <- rep(mean_Y_1, n1)
    # data set 2
    res$alpha_2 <- rep(0, p)
    res$post_mean2 <- rep(0, p)
    res$Xb_2 <- rep(mean_Y_2, n2)
  }
  # return results
  return(res)
}

# #### check results
# time1 <- Sys.time()
# res <- sum_single_effect_two(X_1, Y_1, X_2, Y_2, L = p_c + p_1 + p_2 + 1, scale_x = TRUE, intercept = TRUE)
# Sys.time() - time1
#
# source("two_data_sets/sum_single_effect_two.R")
# time1 <- Sys.time()
# res <- sum_single_effect_two(X_1, Y_1, X_2, Y_2, L = p_c + p_1 + p_2 + 1, scale_x = TRUE, intercept = TRUE)
# Sys.time() - time1