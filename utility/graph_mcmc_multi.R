## MCMC method for Graph
# dta_list are n x p data set
# scale_x : scale the data
# intercept: calculate the mean of Y
# order_int is the initialized order for nodes
# iter_max is the maximun mcmc step
# sigma02_int is initialization for signal prior variance
# sigma2_int is initialization for error variance
# prior_vec : prior for different models
# itermax is the maximum iteration
# tol is the threshold for ELBO
# sigma0_low_bd is the threshold for select effect l

source("utility/graph_given_order_multi.R")
Graph_MCMC_multi <- function(dta_list, scale_x = FALSE, intercept = TRUE, com_mat = NULL,
                             order_int = NULL, iter_max = 50000, sigma02_int = NULL, sigma2_int = NULL,
                             prior_vec = NULL, itermax = 100, L_max = 10, tol = 1e-4, sigma0_low_bd = 1e-8,
                             burn_in = iter_max - 5000, residual_variance_lowerbound = NULL) {
  ## Initialization
  K <- length(dta_list)
  n_group <- 2^K - 1
  p <- rep(NA, K)
  n <- rep(NA, K)
  for (iter_K in seq_len(K)) {
    p[iter_K] <- ncol(dta_list[[iter_K]])
    n[iter_K] <- nrow(dta_list[[iter_K]])
  }
  if (length(unique(p)) > 1) {
    stop("The number of features should be same!")
  } else {
    p <- p[1]
  }
  ## define prior vector
  if (is.null(prior_vec)) {
    prior_vec <- rep(1 / (n_group * p^1.5), n_group)
  }
  lprior_vec <- log(prior_vec)
  # combinatorics matrix
  if (is.null(com_mat)) {
    com_list <- list()
    com_mat <- matrix(c(0, 1), ncol = 1)
    for (iter in 2:K) {
      com_mat_copy <- com_mat
      com_mat <- cbind(1, com_mat)
      com_mat_copy <- cbind(0, com_mat_copy)
      com_mat <- rbind(com_mat_copy, com_mat)
    }
    com_mat <- com_mat[-1, ]
    for (iter_com in seq_len(n_group)) {
      com_list[[iter_com]] <- which(com_mat[iter_com, ] == 1)
    }
  }
  # Initialize order
  if (is.null(order_int)) {
    order_old <- sample(seq_len(p), p)
  } else {
    order_old <- order_int
  }
  # change data
  dta_old_list <- list()
  for (iter_K in seq_len(K)) {
    dta_old_list[[iter_K]] <- dta_list[[iter_K]][, order_old]
  }
  ## Generate the first graph
  res_old <- joint_graph_multi(dta_old_list,
    scale_x = scale_x, intercept = intercept,
    sigma02_int = sigma02_int, sigma2_int = sigma2_int, prior_vec = prior_vec,
    com_mat = com_mat, com_list = com_list, itermax = itermax, L_max = L_max,
    tol = tol, sigma0_low_bd = sigma0_low_bd,
    residual_variance_lowerbound = residual_variance_lowerbound
  )
  ## old results
  alpha_res_old <- res_old$alpha_list
  A_res_old <- res_old$A_list
  sigma2_vec_old <- res_old$sigma2_vec
  llike_mat_old <- res_old$llike_mat
  llike_penalty_old <- res_old$llike_penalty
  llike_old <- sum(llike_mat_old) + sum(llike_penalty_old)
  llike_vec <- rep(NA, iter_max)
  ## save lists
  alpha_list <- list()
  A_list <- list()
  order_list <- list()
  for (iter_K in seq_len(K)) {
    alpha_list[[iter_K]] <- matrix(0, nrow = p, ncol = p)
    A_list[[iter_K]] <- matrix(0, nrow = p, ncol = p)
  }
  ## begin iteration
  for (iter_MCMC in seq_len(iter_max)) {
    if (iter_MCMC %% 1000 == 0) print(iter_MCMC)
    ## Initialize proposal
    dta_pro_list <- list()
    for (iter_K in seq_len(K)) {
      dta_pro_list[[iter_K]] <- dta_old_list[[iter_K]]
    }
    order_pro <- order_old
    # log likelihood
    sigma2_vec_pro <- sigma2_vec_old
    llike_mat_pro <- llike_mat_old
    llike_penalty_pro <- llike_penalty_old
    ## propose the new order
    pos_change <- sample(seq_len(p - 1), 1)
    llike_pro <- llike_old - sum(llike_mat_old[c(pos_change, pos_change + 1), ]) -
      sum(llike_penalty_old[c(pos_change, pos_change + 1)])
    for (iter_K in seq_len(K)) {
      dta_pro_list[[iter_K]][, c(pos_change, pos_change + 1)] <- dta_old_list[[iter_K]][, c(pos_change + 1, pos_change)]
    }
    order_pro[c(pos_change, pos_change + 1)] <- order_old[c(pos_change + 1, pos_change)]
    ## doing variable selection
    if (pos_change == 1) {
      res_pos <- list()
      res_pos$res <- list()
      Y_list <- list()
      for (iter_K in seq_len(K)) {
        Y_list[[iter_K]] <- dta_pro_list[[iter_K]][, 1]
        res_pos$res[[iter_K]] <- list()
        res_pos$res[[iter_K]]$alpha <- rep(0, p)
        res_pos$res[[iter_K]]$post_mean <- rep(0, p)
        if (intercept) {
          res_pos$res[[iter_K]]$Xb <- rep(mean(dta_pro_list[[iter_K]][, 1]), n[iter_K])
        } else {
          res_pos$res[[iter_K]]$Xb <- rep(0, n[iter_K])
        }
      }
      res_pos$sigma2 <- var(unlist(Y_list))
    } else {
      dta_tmp_list <- list()
      for (iter_K in seq_len(K)) {
        dta_tmp_list[[iter_K]] <- list()
        dta_tmp_list[[iter_K]]$X <- dta_pro_list[[iter_K]][, seq_len(pos_change - 1), drop = FALSE]
        dta_tmp_list[[iter_K]]$Y <- dta_pro_list[[iter_K]][, pos_change]
      }
      res_pos <- sum_single_effect_multi(
        dta_list = dta_tmp_list, scale_x = scale_x, intercept = intercept,
        sigma02_int = sigma02_int, sigma2_int = sigma2_vec_old[pos_change + 1], prior_vec = prior_vec,
        L = min(pos_change - 1, L_max), itermax = itermax, tol = tol, sigma0_low_bd = sigma0_low_bd,
        residual_variance_lowerbound = residual_variance_lowerbound
      )
      llike_penalty_pro[pos_change] <- sum(res_pos$alpha * rep(lprior_vec, each = pos_change - 1))
    }
    dta_tmp_list <- list()
    for (iter_K in seq_len(K)) {
      dta_tmp_list[[iter_K]] <- list()
      dta_tmp_list[[iter_K]]$X <- dta_pro_list[[iter_K]][, seq_len(pos_change), drop = FALSE]
      dta_tmp_list[[iter_K]]$Y <- dta_pro_list[[iter_K]][, pos_change + 1]
    }
    res_pos1 <- sum_single_effect_multi(
      dta_list = dta_tmp_list, scale_x = scale_x, intercept = intercept,
      sigma02_int = sigma02_int, sigma2_int = sigma2_vec_old[pos_change], prior_vec = prior_vec,
      L = min(pos_change, L_max), itermax = itermax, tol = tol, sigma0_low_bd = sigma0_low_bd,
      residual_variance_lowerbound = residual_variance_lowerbound
    )
    llike_penalty_pro[pos_change + 1] <- sum(res_pos1$alpha * rep(lprior_vec, each = pos_change))
    # likelihood
    sigma2_vec_pro[c(pos_change, pos_change + 1)] <- c(res_pos$sigma2, res_pos1$sigma2)
    for (iter_K in seq_len(K)) {
      llike_mat_pro[pos_change, iter_K] <- sum(dnorm(
        x = dta_pro_list[[iter_K]][, pos_change],
        mean = res_pos$res[[iter_K]]$Xb, sd = sqrt(res_pos$sigma2), log = TRUE
      ))
      llike_mat_pro[pos_change + 1, iter_K] <- sum(dnorm(
        x = dta_pro_list[[iter_K]][, pos_change + 1],
        mean = res_pos1$res[[iter_K]]$Xb, sd = sqrt(res_pos1$sigma2), log = TRUE
      ))
    }
    llike_pro <- llike_pro + sum(llike_mat_pro[c(pos_change, pos_change + 1), ]) +
      sum(llike_penalty_pro[c(pos_change, pos_change + 1)])

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
      for (iter_K in seq_len(K)) {
        alpha_res_old[[iter_K]][, c(pos_change, pos_change + 1)] <- alpha_res_old[[iter_K]][, c(pos_change + 1, pos_change)]
        A_res_old[[iter_K]][, c(pos_change, pos_change + 1)] <- A_res_old[[iter_K]][, c(pos_change + 1, pos_change)]
        alpha_res_old[[iter_K]][c(pos_change, pos_change + 1), ] <- 0
        A_res_old[[iter_K]][c(pos_change, pos_change + 1), ] <- 0
        alpha_res_old[[iter_K]][pos_change, seq_len(pos_change - 1)] <- res_pos$res[[iter_K]]$alpha
        alpha_res_old[[iter_K]][pos_change + 1, seq_len(pos_change)] <- res_pos1$res[[iter_K]]$alpha
        A_res_old[[iter_K]][pos_change, seq_len(pos_change - 1)] <- res_pos$res[[iter_K]]$post_mean
        A_res_old[[iter_K]][pos_change + 1, seq_len(pos_change)] <- res_pos1$res[[iter_K]]$post_mean
      }
      # likelihood
      sigma2_vec_old <- sigma2_vec_pro
      llike_mat_old <- llike_mat_pro
      llike_penalty_old <- llike_penalty_pro
      llike_old <- llike_pro
      # data and order
      dta_old_list <- dta_pro_list
      order_old <- order_pro
    }
    # save lists
    llike_vec[iter_MCMC] <- llike_old
    if (iter_MCMC > burn_in) {
      order_list[[iter_MCMC - burn_in]] <- order_old
      order_tmp <- order(order_old)
      for (iter_K in seq_len(K)) {
        alpha_list[[iter_K]] <- alpha_list[[iter_K]] +
          alpha_res_old[[iter_K]][order_tmp, order_tmp]
        A_list[[iter_K]] <- A_list[[iter_K]] +
          A_res_old[[iter_K]][order_tmp, order_tmp]
      }
    }
  }

  for (iter_K in seq_len(K)) {
    alpha_list[[iter_K]] <- alpha_list[[iter_K]] / length(order_list)
    A_list[[iter_K]] <- A_list[[iter_K]] / length(order_list)
  }
  # return results
  return(list(
    alpha_list = alpha_list, A_list = A_list,
    llike_vec = llike_vec, order_list = order_list
  ))
}
