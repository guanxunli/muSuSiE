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
iter_max <- 5e5
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