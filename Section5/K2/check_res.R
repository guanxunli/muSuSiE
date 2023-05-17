library(pcalg)
source("Section5/graph_generation.R")
p <- 100
n_tol <- 600
K <- 2
n <- n_tol / K
e_com <- 50 # 50 100 100
e_pri <- 50 # 50 50 20

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

## define metric function
## remove order edge
check_edge <- function(adj_pre, adj_act) {
  adj_pre <- ceiling((adj_pre + t(adj_pre)) / 2)
  adj_act <- ceiling((adj_act + t(adj_act)) / 2)
  return(sum(abs(adj_pre - adj_act)) / 2)
}
## True positive rate
TPrate_fun <- function(adj_pre, adj_act) {
  adj_pre <- ceiling((adj_pre + t(adj_pre)) / 2)
  adj_act <- ceiling((adj_act + t(adj_act)) / 2)
  P <- which(adj_act == 1)
  PP <- which(adj_pre == 1)
  return(length(intersect(P, PP)) / length(P))
}
## False positive rate
FPrate_fun <- function(adj_pre, adj_act) {
  adj_pre <- ceiling((adj_pre + t(adj_pre)) / 2)
  adj_act <- ceiling((adj_act + t(adj_act)) / 2)
  N <- which(adj_act == 0)
  PP <- which(adj_pre == 1)
  return(length(intersect(N, PP)) / length(N))
}
## check adjacency matrix
check_adj_l2 <- function(adj_pre, adj_act) {
  adj_pre <- ceiling((adj_pre + t(adj_pre)) / 2)
  adj_act <- ceiling((adj_act + t(adj_act)) / 2)
  return(sum((adj_pre - adj_act)^2) / 2)
}
check_adj_mcmc <- function(adj_pre, adj_act) {
  adj_pre <- adj_pre + t(adj_pre)
  adj_act <- adj_act + t(adj_act)
  return(sum((adj_pre - adj_act)^2) / 2)
}

############################# PC method ##########################
alphas <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)
out_res <- readRDS(paste0("Section5/K2/results/pc_com", e_com, "pri", e_pri, ".rds"))
## check results
res <- list()
for (iter_K in seq_len(K)) {
  res[[iter_K]] <- list()
}
for (iter_alpha in seq_len(length(alphas))) {
  for (iter_K in seq_len(K)) {
    res[[iter_K]][[iter_alpha]] <- matrix(NA, nrow = n_graph, ncol = 4)
  }
  for (iter_graph in seq_len(n_graph)) {
    for (iter_K in seq_len(K)) {
      ## load true value
      adj_true <- t(graph_sim$G[[iter_graph]][[iter_K]])
      g_true <- as(getGraph(adj_true), "graphNEL")
      ## load results
      adj <- out_res[[iter_graph]][[1]][[iter_K]][[iter_alpha]]
      g <- as(adj, "graphNEL")
      res[[iter_K]][[iter_alpha]][iter_graph, ] <- c(
        check_edge(adj_true, adj),
        TPrate_fun(adj_pre = adj, adj_act = adj_true),
        FPrate_fun(adj_pre = adj, adj_act = adj_true),
        check_adj_l2(adj_pre = adj, adj_act = adj_true)
      )
    }
  }
}

res_ave <- list()
for (iter_alpha in seq_len(length(alphas))) {
  res_ave[[iter_alpha]] <- matrix(0, nrow = n_graph, ncol = 4)
  for (iter_K in seq_len(K)) {
    res_ave[[iter_alpha]] <- res_ave[[iter_alpha]] + res[[iter_K]][[iter_alpha]]
  }
  res_ave[[iter_alpha]] <- res_ave[[iter_alpha]] / K
}

## print results
for (iter_alpha in seq_len(length(alphas))) {
  res_tmp <- round(colMeans(res_ave[[iter_alpha]]), 4)
  cat(
    K, "&", alphas[iter_alpha], "&", e_com, "&", e_pri, "&",
    res_tmp[1], "&", res_tmp[2], "&", res_tmp[3], "&", res_tmp[4], "\\\\", "\n"
  )
}

############################# GES method #############################
lambdas <- c(1, 2, 3, 4, 5)
out_res <- readRDS(paste0("Section5/K2/results/GES_com", e_com, "pri", e_pri, ".rds"))
## check results
res <- list()
for (iter_K in seq_len(K)) {
  res[[iter_K]] <- list()
}
for (iter_lambda in seq_len(length(lambdas))) {
  for (iter_K in seq_len(K)) {
    res[[iter_K]][[iter_lambda]] <- matrix(NA, nrow = n_graph, ncol = 4)
  }
  for (iter_graph in seq_len(n_graph)) {
    for (iter_K in seq_len(K)) {
      ## load true value
      adj_true <- t(graph_sim$G[[iter_graph]][[iter_K]])
      g_true <- as(getGraph(adj_true), "graphNEL")
      ## load results
      adj <- out_res[[iter_graph]][[1]][[iter_K]][[iter_lambda]]
      g <- as(adj, "graphNEL")
      res[[iter_K]][[iter_lambda]][iter_graph, ] <- c(
        check_edge(adj_true, adj),
        TPrate_fun(adj_pre = adj, adj_act = adj_true),
        FPrate_fun(adj_pre = adj, adj_act = adj_true),
        check_adj_l2(adj_pre = adj, adj_act = adj_true)
      )
    }
  }
}

res_ave <- list()
for (iter_lambda in seq_len(length(lambdas))) {
  res_ave[[iter_lambda]] <- matrix(0, nrow = n_graph, ncol = 4)
  for (iter_K in seq_len(K)) {
    res_ave[[iter_lambda]] <- res_ave[[iter_lambda]] + res[[iter_K]][[iter_lambda]]
  }
  res_ave[[iter_lambda]] <- res_ave[[iter_lambda]] / K
}

## print results
for (iter_lambda in seq_len(length(lambdas))) {
  res_tmp <- round(colMeans(res_ave[[iter_lambda]]), 4)
  cat(
    K, "&", lambdas[iter_lambda], "&", e_com, "&", e_pri, "&",
    res_tmp[1], "&", res_tmp[2], "&", res_tmp[3], "&", res_tmp[4], "\\\\", "\n"
  )
}

############################# joint GES method #############################
lambdas <- c(1, 2, 3, 4, 5)
out_res <- readRDS(paste0("Section5/K2/results/jointGES_com", e_com, "pri", e_pri, ".rds"))
## check results
res <- list()
for (iter_K in seq_len(K)) {
  res[[iter_K]] <- list()
}
for (iter_lambda in seq_len(length(lambdas))) {
  for (iter_K in seq_len(K)) {
    res[[iter_K]][[iter_lambda]] <- matrix(NA, nrow = n_graph, ncol = 4)
  }
  for (iter_graph in seq_len(n_graph)) {
    for (iter_K in seq_len(K)) {
      ## load true value
      adj_true <- t(graph_sim$G[[iter_graph]][[iter_K]])
      g_true <- as(getGraph(adj_true), "graphNEL")
      ## load results
      adj <- out_res[[iter_graph]][[1]][[iter_K]][[iter_lambda]]
      g <- as(adj, "graphNEL")
      res[[iter_K]][[iter_lambda]][iter_graph, ] <- c(
        check_edge(adj_true, adj),
        TPrate_fun(adj_pre = adj, adj_act = adj_true),
        FPrate_fun(adj_pre = adj, adj_act = adj_true),
        check_adj_l2(adj_pre = adj, adj_act = adj_true)
      )
    }
  }
}

res_ave <- list()
for (iter_lambda in seq_len(length(lambdas))) {
  res_ave[[iter_lambda]] <- matrix(0, nrow = n_graph, ncol = 4)
  for (iter_K in seq_len(K)) {
    res_ave[[iter_lambda]] <- res_ave[[iter_lambda]] + res[[iter_K]][[iter_lambda]]
  }
  res_ave[[iter_lambda]] <- res_ave[[iter_lambda]] / K
}

## print results
for (iter_lambda in seq_len(length(lambdas))) {
  res_tmp <- round(colMeans(res_ave[[iter_lambda]]), 4)
  cat(
    K, "&", lambdas[iter_lambda], "&", e_com, "&", e_pri, "&",
    res_tmp[1], "&", res_tmp[2], "&", res_tmp[3], "&", res_tmp[4], "\\\\", "\n"
  )
}

############################# muSuSiE-DAG #############################
for (iter_prior in seq_len(4)) {
  out_res <- readRDS(paste0(
    "Section5/K2/results/muSuSiE_prior", iter_prior,
    "com", e_com, "pri", e_pri, ".rds"
  ))
  ## check results
  res <- list()
  for (iter_K in seq_len(K)) {
    res[[iter_K]] <- matrix(NA, nrow = n_graph, ncol = 4)
  }
  res_ave <- matrix(0, nrow = n_graph, ncol = 4)

  for (iter_graph in seq_len(n_graph)) {
    res_tmp <- out_res[[iter_graph]][[1]]
    A_mat_list <- res_tmp$A_list
    alpha_mat_list <- res_tmp$alpha_list
    for (iter_K in seq_len(K)) {
      adj <- ifelse(alpha_mat_list[[iter_K]] > 0.5, 1, 0)
      adj <- t(adj)
      g <- as(getGraph(adj), "graphNEL")
      adj_true <- t(graph_sim$G[[iter_graph]][[iter_K]])
      g_true <- as(getGraph(adj_true), "graphNEL")
      ## save results
      res[[iter_K]][iter_graph, ] <- c(
        check_edge(adj_true, adj),
        TPrate_fun(adj_pre = adj, adj_act = adj_true),
        FPrate_fun(adj_pre = adj, adj_act = adj_true),
        check_adj_mcmc(adj_pre = alpha_mat_list[[iter_K]], adj_act = adj_true)
      )
    }
  }
  # show results
  for (iter_K in seq_len(K)) {
    res_ave <- res_ave + res[[iter_K]]
  }
  # average results
  res_ave <- round(colMeans(res_ave) / K, 4)
  cat(
    K, "&", "prior", iter_prior, "&", e_com, "&", e_pri, "&",
    res_ave[1], "&", res_ave[2], "&", res_ave[3], "&", res_ave[4], "\\\\", "\n"
  )
}

############################# KSEC #############################
out_res <- readRDS( paste0("Section5/K2/results/JESC_com",
                           e_com, "pri", e_pri, ".rds"))
## check results
res <- list()
for (iter_K in seq_len(K)) {
  res[[iter_K]] <- matrix(NA, nrow = n_graph, ncol = 4)
}
res_ave <- matrix(0, nrow = n_graph, ncol = 4)

for (iter_graph in seq_len(n_graph)) {
  alpha_mat_list <- out_res[[iter_graph]][[1]]
  for (iter_K in seq_len(K)) {
    adj <- ifelse(alpha_mat_list[[iter_K]] > 0.5, 1, 0)
    adj <- t(adj)
    g <- as(getGraph(adj), "graphNEL")
    adj_true <- t(graph_sim$G[[iter_graph]][[iter_K]])
    g_true <- as(getGraph(adj_true), "graphNEL")
    ## save results
    res[[iter_K]][iter_graph, ] <- c(
      check_edge(adj_true, adj),
      TPrate_fun(adj_pre = adj, adj_act = adj_true),
      FPrate_fun(adj_pre = adj, adj_act = adj_true),
      check_adj_mcmc(adj_pre = alpha_mat_list[[iter_K]], adj_act = adj_true)
    )
  }
}
# show results
for (iter_K in seq_len(K)) {
  res_ave <- res_ave + res[[iter_K]]
}
# average results
res_ave <- round(colMeans(res_ave) / K, 4)
cat(
  K, "&", e_com, "&", e_pri, "&",
  res_ave[1], "&", res_ave[2], "&", res_ave[3], "&", res_ave[4], "\\\\", "\n"
)
