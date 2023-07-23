library(pcalg)
library(parallel)
source("Section5/graph_generation.R")
p <- 100
n_tol <- 600
K <- 2
n <- n_tol / K
e_com <- 50
e_pri <- 50

####### generate graphs ##############
library(foreach)
library(doParallel)
library(doRNG)
#### generate graph
set.seed(2021)
n_graph <- 50
graph_sim <- graph_generation(
  K = K, n_graph = n_graph, p = p, n_tol = n_tol,
  e_com = e_com, e_pri = e_pri
)

############### muSuSiE ###########
source("Section5/mcmc_con/check_con.R")
#### define prior
prior_pi <- c(1 / (2 * p^2), 1 / p^2.25)

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
prior_vec <- rep(NA, n_group)
for (iter in seq_len(n_group)) {
  prior_vec[iter] <- prior_pi[com_length[[iter]]]
}

#### simulations
iter_max <- 1e6
n_simu <- 50
out_res <- list()
## do parallel
set.seed(2021)
cl <- makeCluster(50)
registerDoParallel(cl)
out_res <- foreach(iter = seq_len(n_simu)) %dorng% {
  library(pcalg)
  # get order
  dta <- matrix(NA, nrow = K * n, ncol = p)
  adj_true <- list()
  for (iter_K in seq_len(K)) {
    dta[(1 + (iter_K - 1) * n):(iter_K * n), ] <- graph_sim$X[[1]][[iter_K]]
    adj_true[[iter_K]] <- t(graph_sim$G[[1]][[iter_K]])
  }
  score_ges <- new("GaussL0penObsScore", data = dta, intercept = FALSE)
  ges_fit <- ges(score_ges)
  ges_adj <- as(ges_fit$repr, "matrix")
  ges_adj <- ifelse(ges_adj == TRUE, 1, 0)
  graph_i <- igraph::graph_from_adjacency_matrix(ges_adj, mode = "directed", diag = FALSE)
  order_int <- as.numeric(igraph::topo_sort(graph_i))
  # Do MCMC
  res <- Graph_MCMC_multi(
    dta_list = graph_sim$X[[1]], scale_x = FALSE, intercept = TRUE,
    order_int = order_int, iter_max = iter_max,
    prior_vec = prior_vec, itermax = 100, L_max = 10
  )
  return(res)
}
stopCluster(cl)
saveRDS(out_res, paste0("Section5/mcmc_con/results/muSuSiE_com", e_com, "pri", e_pri, ".rds"))

#### check convergence
library(ggplot2)
out_res <- readRDS(paste0("Section5/mcmc_con/results/muSuSiE_com", e_com, "pri", e_pri, ".rds"))
loglikelihood_mat <- matrix(as.numeric(unlist(out_res)), nrow = 10001, ncol = n_simu)
sd_loglikelihood <- apply(loglikelihood_mat, 1, sd)
mean_loglikelihood <- apply(loglikelihood_mat, 1, mean)
df_plot <- data.frame(
  loglikelihood = mean_loglikelihood,
  sd_loglikelihood = sd_loglikelihood,
  iteration = seq_len(10001)
)
p1 <- ggplot(df_plot, aes(x = iteration, y = loglikelihood)) +
  geom_line(linewidth = 1.5) +
  geom_line(aes(x = iteration, y = loglikelihood - sd_loglikelihood), linetype = "dashed", linewidth = 0.5) +
  geom_line(aes(x = iteration, y = loglikelihood + sd_loglikelihood), linetype = "dashed", linewidth = 0.5) +
  ylab("Log likelihood") +
  xlab("Number of iterations")
for (iter in seq_len(50)) {
  p1 <- p1  + geom_line(aes(x = iteration), y = loglikelihood_mat[, iter], linewidth = 0.1, alpha = 0.1)
}
p1 <- p1 + scale_x_continuous(
  breaks = c(0, 2e3, 4e3, 6e3, 8e3, 1e4),
  labels = c("0K", "200K", "400K", "600K", "800K", "1000K")
) +
  theme(
    axis.title = element_text(size = 24),
    panel.grid.minor = element_line(linewidth = 0.2),
    panel.grid.major = element_line(linewidth = 0.5),
    legend.position = "none",
    legend.key.size = unit(2.3, "lines"),
    legend.box.margin = margin(t = 0, r = 0),
    legend.title = element_blank()
  ) +
  theme_bw(base_size = 24)
ggsave(paste0("Section5/mcmc_con/results/", e_com, "pri", e_pri, ".pdf"),
       p1,
       width = 10, height = 5
)
