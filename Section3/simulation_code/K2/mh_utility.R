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
mh_joint <- function(X1, Y1, X2, Y2, state_init = NULL, max_mcmc = 1e5, burn_in = 1e4,
                     kappa0 = 2, kappa1 = 1.5) {
  # Initialization
  n <- nrow(X1)
  p <- ncol(X1)
  if (is.null(state_init)) {
    state_init <- rep(0, p)
  }
  state_curr1 <- state_curr2 <- state_init
  sum_Y1 <- sum(Y1^2)
  sum_Y2 <- sum(Y2^2)
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
  logpost_curr1 <- logpost_fun(X1, Y1, sum_Y1, state_curr1)
  logpost_curr2 <- logpost_fun(X2, Y2, sum_Y2, state_curr2)
  logpost_curr <- logpost_curr1 + logpost_curr2
  ## begin iteration
  est_pip1 <- numeric(p)
  est_pip2 <- numeric(p)
  for (iter_mcmc in seq_len(max_mcmc)) {
    ## propose new state
    index_prop <- sample(p, 1)
    dta_prop <- sample(c(1, 2, 3), size = 1, prob = c(1 / 3, 1 / 3, 1 / 3))
    state_prop1 <- state_curr1
    state_prop2 <- state_curr2
    if (dta_prop == 1) {
      state_prop1[index_prop] <- 1 - state_prop1[index_prop]
      logpost_prop1 <- logpost_fun(X1, Y1, sum_Y1, state_prop1)
      logpost_prop2 <- logpost_curr2
    } else if (dta_prop == 2) {
      state_prop2[index_prop] <- 1 - state_prop2[index_prop]
      logpost_prop2 <- logpost_fun(X2, Y2, sum_Y2, state_prop2)
      logpost_prop1 <- logpost_curr1
    } else {
      state_prop1[index_prop] <- 1 - state_prop1[index_prop]
      logpost_prop1 <- logpost_fun(X1, Y1, sum_Y1, state_prop1)
      state_prop2[index_prop] <- 1 - state_prop2[index_prop]
      logpost_prop2 <- logpost_fun(X2, Y2, sum_Y2, state_prop2)
    }
    logpost_prop <- logpost_prop1 + logpost_prop2
    accept_prob <- min(1, (exp(logpost_prop - logpost_curr)))
    ## accept or nor
    if (runif(1, 0, 1) < accept_prob) {
      state_curr1 <- state_prop1
      logpost_curr1 <- logpost_prop1
      state_curr2 <- state_prop2
      logpost_curr2 <- logpost_prop2
    }
    if (iter_mcmc > burn_in) {
      est_pip1 <- est_pip1 + state_curr1
      est_pip2 <- est_pip2 + state_curr2
    }
  }
  est_pip1 <- est_pip1 / (max_mcmc - burn_in)
  est_pip2 <- est_pip2 / (max_mcmc - burn_in)
  return(list(est_pip1 = est_pip1, est_pip2 = est_pip2))
}