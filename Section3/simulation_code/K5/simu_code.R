set.seed(1)
library(foreach)
library(doParallel)
library(doRNG)
#### Initialization
## Define parameters
source("utility/sum_single_effect_multi.R")
## simulation function
simu_vs_fun <- function(K, n, p, p_c, p_s, sigma, sigma0, prior_vec) {
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
  time1 <- Sys.time()
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
  time_c <- Sys.time() - time1
  #### Single data set
  index_res_s <- list()
  sens_res_s <- rep(NA, K)
  prec_res_s <- rep(NA, K)
  ## data set
  time1 <- Sys.time()
  for (k in seq_len(K)) {
    res <- susieR::susie(X = dta_list[[k]]$X, y = dta_list[[k]]$Y, L = p_c + p_s + 1)
    index_res_s[[k]] <- as.numeric(res$sets$cs)
    if (length(index_res_s[[k]]) == 0) {
      sens_res_s[k] <- prec_res_s[k] <- 0
    } else {
      sens_res_s[k] <- length(intersect(
        index_res_s[[k]],
        c(index_s[[k]], index_c)
      )) / (p_c + p_s)
      prec_res_s[k] <- length(intersect(
        index_res_s[[k]],
        c(index_s[[k]], index_c)
      )) / length(index_res_s[[k]])
    }
  }
  time_s <- Sys.time() - time1
  #### LASSO method
  index_res_l <- list()
  sens_res_l <- rep(NA, K)
  prec_res_l <- rep(NA, K)
  ## data set
  time1 <- Sys.time()
  for (k in seq_len(K)) {
    cv_model <- glmnet::cv.glmnet(x = dta_list[[k]]$X, y = dta_list[[k]]$Y, alpha = 1)
    lambda_cv <- cv_model$lambda.min
    res <- glmnet::glmnet(x = dta_list[[k]]$X, y = dta_list[[k]]$Y, alpha = 1, lambda = lambda_cv)
    index_res_l[[k]] <- which(as.numeric(res$beta) != 0)
    if (length(index_res_l[[k]]) == 0) {
      sens_res_l[k] <- prec_res_l[k] <- 0
    } else {
      sens_res_l[k] <- length(intersect(
        index_res_l[[k]],
        c(index_s[[k]], index_c)
      )) / (p_c + p_s)
      prec_res_l[k] <- length(intersect(
        index_res_l[[k]],
        c(index_s[[k]], index_c)
      )) / length(index_res_l[[k]])
    }
  }
  time_l <- Sys.time() - time1
  #### return
  return(list(
    sens_res_c = mean(sens_res_c), sens_res_s = mean(sens_res_s),
    sens_res_l = mean(sens_res_l), prec_res_c = mean(prec_res_c),
    prec_res_s = mean(prec_res_s), prec_res_l = mean(prec_res_l),
    time_c = time_c, time_s = time_s, time_l = time_l
  ))
}

############################## do simulations ###############################
K <- 5
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
set_df <- data.frame(n = c(100, 500), p = c(600, 1000))
sigma_vec <- c(1, 2)
p_c_vec <- c(10, 25)
p_s_vec <- c(2, 5)
iter_sim_max <- 500

## do simulations
for (iter_set in seq_len(nrow(set_df))) {
  n <- set_df$n[iter_set]
  p <- set_df$p[iter_set]
  prior_vec <- rep(NA, n_group)
  prior_pi <- c(1 / p^1.4, 1 / p^1.55, 1 / p^1.7, 1 / p^1.85, 1 / p^2)
  for (iter in seq_len(n_group)) {
    prior_vec[iter] <- prior_pi[com_length[[iter]]]
  }
  for (sigma in sigma_vec) {
    for (p_c in p_c_vec) {
      for (p_s in p_s_vec) {
        cl <- makeCluster(25)
        registerDoParallel(cl)
        set.seed(2022)
        out_res <- foreach(iter = seq_len(iter_sim_max)) %dorng% {
          simu_vs_fun(
            K = K, n = n, p = p, p_c = p_c, p_s = p_s,
            sigma = sigma, sigma0 = sigma0, prior_vec = prior_vec
          )
        }
        stopCluster(cl)
        print(paste0(
          "Finish K", K, "n", n, "p", p, "sigma", sigma,
          "pc", p_c, "ps", p_s
        ))
        saveRDS(out_res, paste0(
          "Section3/simulation_code/K5/results/K", K,
          "n", n, "p", p, "sigma", sigma, "pc", p_c,
          "ps", p_s, ".rds"
        ))
      }
    }
  }
}

############################### check results ###############################
library(ggplot2)
library(gridExtra)
## Compare muSuSiE and SuSiE
for (iter_set in seq_len(nrow(set_df))) {
  n <- set_df$n[iter_set]
  p <- set_df$p[iter_set]
  for (sigma in sigma_vec) {
    for (p_c in p_c_vec) {
      for (p_s in p_s_vec) {
        out_res <- readRDS(paste0(
          "Section3/simulation_code/K2/results/K", K, "n", n,
          "p", p, "sigma", sigma, "pc", p_c, "ps", p_s, ".rds"
        ))
        n_iter <- length(out_res)
        out_df <- matrix(unlist(out_res),
          nrow = n_iter, ncol = length(out_res[[1]]),
          byrow = TRUE
        )
        cat(
          p, "&", n, "&", p_c, "&", p_s, "&", sigma^2, "&",
          round(colMeans(out_df), 4)[c(1, 2, 4, 5)], "\n"
        )
      }
    }
  }
}

## results for LASSO
for (iter_set in seq_len(nrow(set_df))) {
  n <- set_df$n[iter_set]
  p <- set_df$p[iter_set]
  for (sigma in sigma_vec) {
    for (p_c in p_c_vec) {
      for (p_s in p_s_vec) {
        out_res <- readRDS(paste0(
          "Section3/simulation_code/K2/results/K", K, "n", n,
          "p", p, "sigma", sigma, "pc", p_c, "ps", p_s, ".rds"
        ))
        n_iter <- length(out_res)
        out_df <- matrix(unlist(out_res),
          nrow = n_iter, ncol = length(out_res[[1]]),
          byrow = TRUE
        )
        cat(
          p, "&", n, "&", p_c, "&", p_s, "&", sigma^2, "&",
          round(colMeans(out_df), 4)[c(3, 6)], "\n"
        )
      }
    }
  }
}

## compare time
for (iter_set in seq_len(nrow(set_df))) {
  n <- set_df$n[iter_set]
  p <- set_df$p[iter_set]
  for (sigma in sigma_vec) {
    for (p_c in p_c_vec) {
      for (p_s in p_s_vec) {
        out_res <- readRDS(paste0(
          "Section3/simulation_code/K2/results/K", K, "n", n,
          "p", p, "sigma", sigma, "pc", p_c, "ps", p_s, ".rds"
        ))
        n_iter <- length(out_res)
        out_df <- matrix(unlist(out_res),
          nrow = n_iter, ncol = length(out_res[[1]]),
          byrow = TRUE
        )
        cat(
          p, "&", n, "&", p_c, "&", p_s, "&", sigma^2, "&",
          round(colMeans(out_df), 4)[c(7, 8)], "\n"
        )
      }
    }
  }
}

############################### plot one figure ###############################
## define parameters
K <- 2
sigma <- 1
p_c <- 10
p_s <- 2
# first two figures
n <- 100
p <- 600
out_res <- readRDS(paste0(
  "Section3/simulation_code/K2/results/K", K, "n", n,
  "p", p, "sigma", sigma, "pc", p_c, "ps", p_s, ".rds"
))
n_iter <- length(out_res)
out_df <- matrix(unlist(out_res),
  nrow = n_iter, ncol = length(out_res[[1]]),
  byrow = TRUE
)[, c(1, 2, 4, 5)]
colnames(out_df) <- c("sens_mu", "sens_si", "prec_mu", "prec_si")
out_df <- as.data.frame(out_df)

## plot sensitivity
out_sens <- data.frame(
  sens = c(out_df$sens_mu, out_df$sens_si),
  group = c(rep("multiple", n_iter), rep("single", n_iter))
)
p1 <- ggplot(out_sens, aes(x = group, y = sens, fill = group)) +
  geom_boxplot() +
  labs(title = "sensitivity for p = 600 and  n = 100") +
  xlab("method") +
  ylab("sensitivity") +
  ylim(c(0, 1)) +
  scale_colour_manual(values = c("red", "green4"), breaks = c("multiple", "single")) +
  theme_bw(base_size = 22) +
  theme(legend.position = "none")
# plot precision
out_prec <- data.frame(
  prec = c(out_df$prec_mu, out_df$prec_si),
  group = c(rep("multiple", n_iter), rep("single", n_iter))
)
p2 <- ggplot(out_prec, aes(x = group, y = prec, fill = group)) +
  geom_boxplot() +
  labs(title = "precision for p = 600 and  n = 100") +
  xlab("method") +
  ylab("precision") +
  ylim(c(0, 1)) +
  scale_colour_manual(values = c("red", "green4"), breaks = c("multiple", "single")) +
  theme_bw(base_size = 22) +
  theme(legend.position = "none")

# second two figures
n <- 500
p <- 1000
out_res <- readRDS(paste0(
  "Section3/simulation_code/K2/results/K", K, "n", n,
  "p", p, "sigma", sigma, "pc", p_c, "ps", p_s, ".rds"
))
n_iter <- length(out_res)
out_df <- matrix(unlist(out_res),
  nrow = n_iter, ncol = length(out_res[[1]]),
  byrow = TRUE
)[, c(1, 2, 4, 5)]
colnames(out_df) <- c("sens_mu", "sens_si", "prec_mu", "prec_si")
out_df <- as.data.frame(out_df)

## plot sensitivity
out_sens <- data.frame(
  sens = c(out_df$sens_mu, out_df$sens_si),
  group = c(rep("multiple", n_iter), rep("single", n_iter))
)
p3 <- ggplot(out_sens, aes(x = group, y = sens, fill = group)) +
  geom_boxplot() +
  labs(title = "sensitivity for p = 1000 and  n = 500") +
  xlab("method") +
  ylab("sensitivity") +
  ylim(c(0, 1)) +
  scale_colour_manual(values = c("red", "green4"), breaks = c("multiple", "single")) +
  theme_bw(base_size = 22) +
  theme(legend.position = "none")
# plot precision
out_prec <- data.frame(
  prec = c(out_df$prec_mu, out_df$prec_si),
  group = c(rep("multiple", n_iter), rep("single", n_iter))
)
p4 <- ggplot(out_prec, aes(x = group, y = prec, fill = group)) +
  geom_boxplot() +
  labs(title = "precision for p = 1000 and  n = 500") +
  xlab("method") +
  ylab("precision") +
  ylim(c(0, 1)) +
  scale_colour_manual(values = c("red", "green4"), breaks = c("multiple", "single")) +
  theme_bw(base_size = 22) +
  theme(legend.position = "none")

p4_legends <- ggplot(out_prec, aes(x = group, y = prec, fill = group)) +
  geom_boxplot() +
  labs(title = "precision for p = 1000 and  n = 500") +
  xlab("method") +
  ylab("precision") +
  ylim(c(0, 1)) +
  scale_colour_manual(values = c("red", "green4"), breaks = c("multiple", "single")) +
  theme_bw(base_size = 22) +
  theme(legend.position = "bottom")
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}
shared_legend <- extract_legend(p4_legends)

## output figure
pdf(file = paste0(
  "Section3/simulation_code/K2/results/K", K, "n", n,
  "p", p, "sigma", sigma, ".pdf"
), width = 16, height = 16)
print(grid.arrange(arrangeGrob(p1, p3, p2, p4, ncol = 2, nrow = 2),
  shared_legend,
  nrow = 2, heights = c(10, 1), widths = 11
))
dev.off()
