source("utility/utility_K2/Graph_given_order_two.R")
Graph_MCMC_two <- function(dta_1, dta_2, scale_x = FALSE, intercept = TRUE,
                           order_int = NULL, iter_max = 50000,
                           sigma02_int = NULL, sigma2_int = NULL, prior_vec = NULL,
                           itermax = 100, L_max = 10, tol = 1e-4, sigma0_low_bd = 1e-8,
                           burn_in = iter_max - 5000, residual_variance_lowerbound = NULL) {
  ## Initialization
  p <- ncol(dta_1)
  if (p != ncol(dta_2)) stop("The number of features should be same!")
  n1 <- nrow(dta_1)
  n2 <- nrow(dta_2)
  ## define prior vector
  if (is.null(prior_vec)) {
    prior_vec <- c(1 / (2 * p^1.5), 1 / p^2)
  }
  lprior_vec <- log(prior_vec)
  # Initialize order
  if (is.null(order_int)) {
    order_old <- sample(seq_len(p), p)
  } else {
    order_old <- order_int
  }
  # change data
  dta_1_old <- dta_1[, order_old]
  dta_2_old <- dta_2[, order_old]
  ## load the main function
  res_old <- joint_graph_fun_two(
    dta_1 = dta_1_old, dta_2 = dta_2_old, scale_x = scale_x, intercept = intercept,
    sigma02_int = sigma02_int, sigma2_int = sigma2_int,
    prior_vec = prior_vec, itermax = itermax,
    L_max = L_max, tol = tol, sigma0_low_bd = sigma0_low_bd,
    residual_variance_lowerbound = residual_variance_lowerbound
  )
  # variable selection
  alpha_res_1_old <- res_old$alpha_res_1
  alpha_res_2_old <- res_old$alpha_res_2
  # posterior of parameters
  A_res_1_old <- res_old$A_res_1
  A_res_2_old <- res_old$A_res_2
  # likelihood
  sigma2_vec_old <- res_old$sigma2_vec
  llike_1_vec_old <- res_old$llike_1_vec
  llike_2_vec_old <- res_old$llike_2_vec
  llike_penalty_vec_old <- res_old$llike_penalty_vec
  llike_old <- sum(llike_1_vec_old) + sum(llike_2_vec_old) +
    sum(llike_penalty_vec_old)
  llike_vec <- rep(NA, iter_max)
  ## save lists
  alpha_mat_1 <- matrix(0, nrow = p, ncol = p)
  alpha_mat_2 <- matrix(0, nrow = p, ncol = p)
  A_mat_1 <- matrix(0, nrow = p, ncol = p)
  A_mat_2 <- matrix(0, nrow = p, ncol = p)
  order_list <- list()
  ## load the function
  for (iter_MCMC in seq_len(iter_max)) {
    if (iter_MCMC %% 1000 == 0) print(iter_MCMC)
    ## Initialize proposal
    dta_1_pro <- dta_1_old
    dta_2_pro <- dta_2_old
    order_pro <- order_old
    # log likelihood
    sigma2_vec_pro <- sigma2_vec_old
    llike_1_vec_pro <- llike_1_vec_old
    llike_2_vec_pro <- llike_2_vec_old
    llike_penalty_vec_pro <- llike_penalty_vec_old
    ## propose the new order
    pos_change <- sample(seq_len(p - 1), 1)
    llike_pro <- llike_old - sum(llike_1_vec_old[c(pos_change, pos_change + 1)]) -
      sum(llike_2_vec_old[c(pos_change, pos_change + 1)]) -
      sum(llike_penalty_vec_old[c(pos_change, pos_change + 1)])
    dta_1_pro[, c(pos_change, pos_change + 1)] <- dta_1_old[, c(pos_change + 1, pos_change)]
    dta_2_pro[, c(pos_change, pos_change + 1)] <- dta_2_old[, c(pos_change + 1, pos_change)]
    order_pro[c(pos_change, pos_change + 1)] <- order_old[c(pos_change + 1, pos_change)]
    ## doing variable selection
    if (pos_change == 1) {
      res_pos <- list()
      res_pos$sigma2 <- var(c(dta_1_pro[, 1], dta_2_pro[, 1]))
      res_pos$alpha_1 <- rep(0, p)
      res_pos$post_mean1 <- rep(0, p)
      res_pos$alpha_2 <- rep(0, p)
      res_pos$post_mean2 <- rep(0, p)
      if (intercept) {
        res_pos$Xb_1 <- rep(mean(dta_1_pro[, 1]), n1)
        res_pos$Xb_2 <- rep(mean(dta_2_pro[, 1]), n2)
      } else {
        res_pos$Xb_1 <- rep(0, n1)
        res_pos$Xb_2 <- rep(0, n2)
      }
    } else {
      res_pos <- sum_single_effect_two(
        X_1 = dta_1_pro[, seq_len(pos_change - 1), drop = FALSE], Y_1 = dta_1_pro[, pos_change],
        X_2 = dta_2_pro[, seq_len(pos_change - 1), drop = FALSE], Y_2 = dta_2_pro[, pos_change],
        scale_x = scale_x, intercept = intercept,
        sigma02_int = sigma02_int, sigma2_int = sigma2_vec_old[pos_change + 1],
        prior_vec = prior_vec, L = min(pos_change - 1, L_max),
        itermax = itermax, tol = tol, sigma0_low_bd = sigma0_low_bd,
        residual_variance_lowerbound = residual_variance_lowerbound
      )
      llike_penalty_vec_pro[pos_change] <- sum(res_pos$alpha * c(rep(lprior_vec[1], 2 * (pos_change - 1)), rep(lprior_vec[2], (pos_change - 1))))
    }
    res_pos1 <- sum_single_effect_two(
      X_1 = dta_1_pro[, seq_len(pos_change), drop = FALSE], Y_1 = dta_1_pro[, pos_change + 1],
      X_2 = dta_2_pro[, seq_len(pos_change), drop = FALSE], Y_2 = dta_2_pro[, pos_change + 1],
      scale_x = scale_x, intercept = intercept,
      sigma02_int = sigma02_int, sigma2_int = sigma2_vec_old[pos_change],
      prior_vec = prior_vec, L = min(pos_change, L_max),
      itermax = itermax, tol = tol, sigma0_low_bd = sigma0_low_bd,
      residual_variance_lowerbound = residual_variance_lowerbound
    )
    llike_penalty_vec_pro[pos_change + 1] <- sum(res_pos1$alpha * c(rep(lprior_vec[1], 2 * pos_change), rep(lprior_vec[2], pos_change)))
    # likelihood
    sigma2_vec_pro[c(pos_change, pos_change + 1)] <- c(res_pos$sigma2, res_pos1$sigma2)
    llike_1_vec_pro[pos_change] <- sum(dnorm(x = dta_1_pro[, pos_change], mean = res_pos$Xb_1, sd = sqrt(res_pos$sigma2), log = TRUE))
    llike_2_vec_pro[pos_change] <- sum(dnorm(x = dta_2_pro[, pos_change], mean = res_pos$Xb_2, sd = sqrt(res_pos$sigma2), log = TRUE))
    llike_1_vec_pro[pos_change + 1] <- sum(dnorm(x = dta_1_pro[, pos_change + 1], mean = res_pos1$Xb_1, sd = sqrt(res_pos1$sigma2), log = TRUE))
    llike_2_vec_pro[pos_change + 1] <- sum(dnorm(x = dta_2_pro[, pos_change + 1], mean = res_pos1$Xb_2, sd = sqrt(res_pos1$sigma2), log = TRUE))
    llike_pro <- llike_pro + sum(llike_1_vec_pro[c(pos_change, pos_change + 1)]) +
      sum(llike_2_vec_pro[c(pos_change, pos_change + 1)]) +
      sum(llike_penalty_vec_pro[c(pos_change, pos_change + 1)])
    # accept or not
    if (llike_pro > llike_old) {
      accept <- TRUE
    } else {
      U <- runif(1)
      thres <- exp(llike_pro - llike_old)
      if (U < thres) {
        accept <- TRUE
      } else {
        accept <- FALSE
      }
    }
    # update
    if (accept) {
      # change matrix order
      alpha_res_1_old[, c(pos_change, pos_change + 1)] <- alpha_res_1_old[, c(pos_change + 1, pos_change)]
      alpha_res_2_old[, c(pos_change, pos_change + 1)] <- alpha_res_2_old[, c(pos_change + 1, pos_change)]
      A_res_1_old[, c(pos_change, pos_change + 1)] <- A_res_1_old[, c(pos_change + 1, pos_change)]
      A_res_2_old[, c(pos_change, pos_change + 1)] <- A_res_2_old[, c(pos_change + 1, pos_change)]
      # proposed matrix
      alpha_res_1_old[c(pos_change, pos_change + 1), ] <- 0
      alpha_res_2_old[c(pos_change, pos_change + 1), ] <- 0
      A_res_1_old[c(pos_change, pos_change + 1), ] <- 0
      A_res_2_old[c(pos_change, pos_change + 1), ] <- 0
      # save matrix
      alpha_res_1_old[pos_change, seq_len(pos_change - 1)] <- res_pos$alpha_1
      alpha_res_2_old[pos_change, seq_len(pos_change - 1)] <- res_pos$alpha_2
      alpha_res_1_old[pos_change + 1, seq_len(pos_change)] <- res_pos1$alpha_1
      alpha_res_2_old[pos_change + 1, seq_len(pos_change)] <- res_pos1$alpha_2
      A_res_1_old[pos_change, seq_len(pos_change - 1)] <- res_pos$post_mean1
      A_res_2_old[pos_change, seq_len(pos_change - 1)] <- res_pos$post_mean2
      A_res_1_old[pos_change + 1, seq_len(pos_change)] <- res_pos1$post_mean1
      A_res_2_old[pos_change + 1, seq_len(pos_change)] <- res_pos1$post_mean2
      # likelihood
      sigma2_vec_old <- sigma2_vec_pro
      llike_1_vec_old <- llike_1_vec_pro
      llike_2_vec_old <- llike_2_vec_pro
      llike_penalty_vec_old <- llike_penalty_vec_pro
      llike_old <- llike_pro
      # data and order
      dta_1_old <- dta_1_pro
      dta_2_old <- dta_2_pro
      order_old <- order_pro
    }
    ## save lists
    llike_vec[iter_MCMC] <- llike_old
    if (iter_MCMC > burn_in) {
      order_list[[iter_MCMC - burn_in]] <- order_old
      order_tmp <- order(order_old)
      alpha_mat_1 <- alpha_mat_1 + alpha_res_1_old[order_tmp, order_tmp]
      alpha_mat_2 <- alpha_mat_2 + alpha_res_2_old[order_tmp, order_tmp]
      A_mat_1 <- A_mat_1 + A_res_1_old[order_tmp, order_tmp]
      A_mat_2 <- A_mat_2 + A_res_2_old[order_tmp, order_tmp]
    }
  }
  alpha_mat_1 <- alpha_mat_1 / length(order_list)
  alpha_mat_2 <- alpha_mat_2 / length(order_list)
  A_mat_1 <- A_mat_1 / length(order_list)
  A_mat_2 <- A_mat_2 / length(order_list)
  # return results
  return(list(
    alpha_mat_1 = alpha_mat_1, alpha_mat_2 = alpha_mat_2,
    A_mat_1 = A_mat_1, A_mat_2 = A_mat_2,
    order_list = order_list, llike_vec = llike_vec
  ))
}
