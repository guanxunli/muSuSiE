#################### Metropolis Hasting seperate ####################
mh_seperate <- function(X, Y, state_init = NULL, max_mcmc = 1e5, burn_in = 1e4,
                        kappa0 = 2, kappa1 = 1.5) {
  # Initialization
  n <- nrow(X)
  p <- ncol(X)
  if (is.null(state_init)) {
    state_init <- rep(0, p)
  }
  state_curr <- state_init
  sum_Y <- sum(Y^2)
  g <- p^(2 * kappa1) - 1
  # log posterior
  logpost_fun <- function(gamma_x) {
    # get parameters
    gamma_l1 <- sum(gamma_x)
    tmp1 <- -(kappa0 + kappa1) * gamma_l1 * log(p)
    # calculate posterior
    if (gamma_l1 > 0) {
      X_tmp <- X[, which(gamma_x == 1), drop = FALSE]
      XtY <- crossprod(X_tmp, Y)
      R2 <- as.numeric(crossprod(XtY, solve(crossprod(X_tmp), XtY))) / sum_Y
      tmp2 <- -n / 2 * log(1 + g * (1 - R2))
    } else {
      tmp2 <- -kappa1 * n * log(p)
    }
    return(tmp1 + tmp2 + kappa1 * n * log(p))
  }

  logpost_curr <- logpost_fun(state_curr)
  ## begin iteration
  est_pip <- numeric(p)
  for (iter_mcmc in seq_len(max_mcmc)) {
    ## propose new state
    index_prop <- sample(p, 1)
    state_prop <- state_curr
    state_prop[index_prop] <- 1 - state_prop[index_prop]
    logpost_prop <- logpost_fun(state_prop)
    accept_prob <- min(1, (exp(logpost_prop - logpost_curr)))
    ## accept or nor
    if (runif(1, 0, 1) < accept_prob) {
      state_curr <- state_prop
      logpost_curr <- logpost_prop
    }
    if (iter_mcmc > burn_in) {
      est_pip <- est_pip + state_curr
    }
  }
  est_pip <- est_pip / (max_mcmc - burn_in)
  return(est_pip)
}

#################### Metropolis Hasting joint ####################
mh_joint <- function(dta_list, max_mcmc = 1e5, burn_in = 1e4,
                     kappa0 = 2, kappa1 = 1.5) {
  # Initialization
  K <- length(dta_list)
  n_group <- 2^K - 1
  n <- nrow(dta_list[[1]]$X)
  p <- ncol(dta_list[[1]]$X)
  state_init <- rep(0, p)
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
  state_curr_list <- list()
  for (iter_K in seq_len(K)) {
    state_curr_list[[iter_K]] <- state_init
  }
  sum_list <- list()
  for (iter_K in seq_len(K)) {
    sum_list[[iter_K]] <- sum(dta_list[[iter_K]]$Y^2)
  }
  g <- p^(2 * kappa1) - 1
  # log posterior
  logpost_fun <- function(X, Y, sum_Y, gamma_x) {
    # get parameters
    gamma_l1 <- sum(gamma_x)
    tmp1 <- -(kappa0 + kappa1) * gamma_l1 * log(p)
    # calculate posterior
    if (gamma_l1 > 0) {
      X_tmp <- X[, which(gamma_x == 1), drop = FALSE]
      XtY <- crossprod(X_tmp, Y)
      R2 <- as.numeric(crossprod(XtY, solve(crossprod(X_tmp), XtY))) / sum_Y
      tmp2 <- -n / 2 * log(1 + g * (1 - R2))
    } else {
      tmp2 <- -kappa1 * n * log(p)
    }
    return(tmp1 + tmp2 + kappa1 * n * log(p))
  }
  logpost_curr_vec <- numeric(K)
  for (iter_K in seq_len(K)) {
    logpost_curr_vec[[iter_K]] <- logpost_fun(
      dta_list[[iter_K]]$X,
      dta_list[[iter_K]]$Y,
      sum_list[[iter_K]],
      state_curr_list[[iter_K]]
    )
  }
  logpost_curr <- sum(logpost_curr_vec)
  ## begin iteration
  est_pip <- list()
  for (iter_K in seq_len(K)) {
    est_pip[[iter_K]] <- numeric(p)
  }
  for (iter_mcmc in seq_len(max_mcmc)) {
    ## propose new state
    state_prop_list <- state_curr_list
    logpost_prop_vec <- logpost_curr_vec
    index_prop <- sample(p, 1)
    gamma_prop <- sample(n_group, 1)
    dta_prop <- com_list[[gamma_prop]]
    for (iter_dta in dta_prop) {
      state_prop_list[[iter_dta]][index_prop] <- 1 - state_prop_list[[iter_dta]][index_prop]
      logpost_prop_vec[iter_dta] <- logpost_fun(
        dta_list[[iter_dta]]$X,
        dta_list[[iter_dta]]$Y,
        sum_list[[iter_dta]],
        state_prop_list[[iter_dta]]
      )
    }
    logpost_prop <- sum(logpost_prop_vec)
    accept_prob <- min(1, (exp(logpost_prop - logpost_curr)))
    ## accept or nor
    if (runif(1, 0, 1) < accept_prob) {
      state_curr_list <- state_prop_list
      logpost_curr_vec <- logpost_prop_vec
    }
    if (iter_mcmc > burn_in) {
      for (iter_K in seq_len(K)) {
        est_pip[[iter_K]] <- est_pip[[iter_K]] + state_curr_list[[iter_K]]
      }
    }
  }
  for (iter_K in seq_len(K)) {
    est_pip[[iter_K]] <- est_pip[[iter_K]] / (max_mcmc - burn_in)
  }
  return(est_pip)
}
