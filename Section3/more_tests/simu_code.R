set.seed(1)
library(foreach)
library(doParallel)
library(doRNG)
#### Initialization
## Define parameters
source("utility/sum_single_effect_multi.R")
## simulation function
simu_vs_fun <- function(K, n, p, p_c, p_s, sigma, sigma0, prior_vec_list) {
  ## Generate data
  index_c <- sample(seq_len(p), size = p_c, replace = FALSE)
  b <- list()
  index_s <- list()
  dta_list <- list()
  for (k in seq_len(K)) {
    dta_list[[k]] <- list()
    index_s[[k]] <- sample(setdiff(seq_len(p), index_c), size = p_s, replace = FALSE)
    b[[k]] <- rep(0, p)
    b[[k]][c(index_c, index_s[[k]])] <- rnorm(p_c + p_s, mean = 0, sd = sigma0)
    dta_list[[k]]$X <- scale(matrix(rnorm(p * n), nrow = n, ncol = p))
    y_tmp <- dta_list[[k]]$X %*% b[[k]]
    dta_list[[k]]$Y <- y_tmp + rnorm(n, sd = sigma)
  }
  #### multiple data sets at the same times
  n_prior_vec <- length(prior_vec_list)
  sens_res_list <- list()
  prec_res_list <- list()
  for (iter_prior in seq_len(n_prior_vec)) {
    prior_vec <- prior_vec_list[[iter_prior]]
    res_c <- sum_single_effect_multi(dta_list,
      scale_x = TRUE, intercept = TRUE,
      sigma02_int = NULL, sigma2_int = NULL,
      L = p_c + K * p_s + K,
      prior_vec = prior_vec
    )
    index_res_c <- list()
    sens_res_c <- rep(NA, K)
    prec_res_c <- rep(NA, K)
    for (k in seq_len(K)) {
      index_res_c[[k]] <- which(res_c$res[[k]]$alpha > 0.5)
      if (length(index_res_c[[k]]) == 0) {
        sens_res_c[k] <- prec_res_c[k] <- 0
      } else {
        sens_res_c[k] <- length(intersect(
          index_res_c[[k]],
          c(index_s[[k]], index_c)
        )) / (p_c + p_s)
        prec_res_c[k] <- length(intersect(
          index_res_c[[k]],
          c(index_s[[k]], index_c)
        )) / length(index_res_c[[k]])
      }
    }
    sens_res_list[[iter_prior]] <- mean(sens_res_c)
    prec_res_list[[iter_prior]] <- mean(prec_res_c)
  }
  return(list(sens_res_list = sens_res_list, prec_res_list = prec_res_list))
}

############################## do simulations ###############################
K <- 2
n_group <- 2^K - 1
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
com_length <- lapply(com_list, length)
## define parameters
sigma0 <- 0.6
n <- 500
p <- 1000
sigma <- 1
p_c <- 25
p_s <- 5
iter_sim_max <- 500
## define different prior
prior_vec_list <- list()
prior_vec_list[[1]] <- c(1 / (2 * p^1.1), 1 / (2 * p^1.1), 1 / p^1.25)
prior_vec_list[[2]] <- c(1 / (2 * p^1.1), 1 / (2 * p^1.1), 1 / p^1.5)
prior_vec_list[[3]] <- c(1 / (p^1.25), 1 / (p^1.25), 1 / p^1.5)
prior_vec_list[[4]] <- c(1 / (p^1.25), 1 / (p^1.25), 1 / p^1.75)
prior_vec_list[[5]] <- c(1 / (p^1.25), 1 / (p^1.25), 1 / p^2)
prior_vec_list[[6]] <- c(1 / (2 * p^1.25), 1 / (2 * p^1.25), 1 / p^1.5)

## do simulations
cl <- makeCluster(25)
registerDoParallel(cl)
set.seed(2022)
out_res <- foreach(iter = seq_len(iter_sim_max)) %dorng% {
  simu_vs_fun(
    K = K, n = n, p = p, p_c = p_c, p_s = p_s,
    sigma = sigma, sigma0 = sigma0, prior_vec_list = prior_vec_list
  )
}
stopCluster(cl)
saveRDS(out_res, paste0("Section3/more_tests/results/out_res.rds"))

############################### check results ###############################
n_prior <- length(prior_vec_list)
library(ggplot2)
library(gridExtra)
out_res <- readRDS(paste0("Section3/more_tests/results/out_res.rds"))
n_iter <- length(out_res)
out_df <- matrix(unlist(out_res),
  nrow = n_iter, ncol = 2 * n_prior,
  byrow = TRUE
)
out_df <- as.data.frame(out_df)

## plot sensitivity
out_sens <- data.frame(
  sens = unlist((out_df[, seq_len(n_prior)])),
  group = rep(paste0("prior", seq_len(n_prior)), each = n_iter)
)
out_sens$group <- as.factor(out_sens$group)
p1 <- ggplot(out_sens, aes(x = group, y = sens, fill = group)) +
  geom_boxplot() +
  labs(title = "Sensitivity for different priors") +
  xlab("") +
  ylab("Sensitivity") +
  theme_bw(base_size = 22) +
  theme(legend.position = "bottom")
ggsave("Section3/more_tests/results/sens_res.pdf", p1, width = 10, height = 8)

## plot precision
out_prec <- data.frame(
  prec = unlist((out_df[, (n_prior + 1):(2 * n_prior)])),
  group = rep(paste0("prior", seq_len(n_prior)), each = n_iter)
)
out_prec$group <- as.factor(out_prec$group)
p2 <- ggplot(out_prec, aes(x = group, y = prec, fill = group)) +
  geom_boxplot() +
  labs(title = "Precision for different priors") +
  xlab("") +
  ylab("Precision") +
  theme_bw(base_size = 22) +
  theme(legend.position = "bottom")
ggsave("Section3/more_tests/results/prec_res.pdf", p2, width = 10, height = 8)
