# # load data
load("Section6/ovarian.rda")
p <- ncol(data[[1]])
library(pcalg)
library(graph)
library(parallel)
library(stabs)

################################ PC method ########################
set.seed(1)
alphas <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)

pc_fun <- function(dta, alpha) {
  p <- ncol(dta)
  dta_cor <- cor(dta)
  alpha <- alphas[iter_alpha]
  pc_fit <- pc(
    suffStat = list(C = dta_cor, n = dim(dta)[1]),
    indepTest = gaussCItest, alpha = alpha,
    labels = sapply(1:p, toString)
  )
  dag <- as(pc_fit@graph, "matrix")
  return(ifelse(dag == TRUE, 1, 0))
}

for (iter_alpha in seq_len(length(alphas))) {
  alpha_use <- alphas[iter_alpha]
  dag1 <- pc_fun(data[[1]], alpha_use)
  dag2 <- pc_fun(data[[2]], alpha_use)
  pc_adj1 <- dag1 | t(dag1)
  n1 <- sum(pc_adj1) / 2
  pc_adj2 <- dag2 | t(dag2)
  n2 <- sum(pc_adj2) / 2
  pc_adj <- pc_adj1 & pc_adj2
  n_com <- sum(pc_adj) / 2
  n_total <- n1 + n2 - n_com
  n_ratio <- n_com / n_total
  ## check results
  cat(
    "PC & $\\alpha = ", alpha_use, "$&", n1, "&", n2, "&", n_com,
    "&", n_total, "&", round(n_ratio, 4), "\\\\\n"
  )
  # cat("alpha: ", alpha_use, c(sum(pc_adj1), sum(pc_adj2), sum(pc_adj)) / 2, "\n")
}

######################### GES method ########################
set.seed(1)
intercept_use <- TRUE
lambdas <- c(1, 2, 3, 4, 5)

ges_fun <- function(dta, lambda) {
  p <- ncol(dta)
  l0score <- new("GaussL0penObsScore", data = dta, lambda = lambda * log(p), intercept = intercept_use)
  ges_fit <- ges(l0score)
  dag <- as(ges_fit$repr, "matrix")
  return(ifelse(dag == TRUE, 1, 0))
}

for (iter_lambda in seq_len(length(lambdas))) {
  lambda_use <- lambdas[iter_lambda]
  dag1 <- ges_fun(data[[1]], lambda_use)
  dag2 <- ges_fun(data[[2]], lambda_use)
  ges_adj1 <- dag1 | t(dag1)
  n1 <- sum(ges_adj1) / 2
  ges_adj2 <- dag2 | t(dag2)
  n2 <- sum(ges_adj2) / 2
  ges_adj <- ges_adj1 & ges_adj2
  n_com <- sum(ges_adj) / 2
  n_total <- n1 + n2 - n_com
  n_ratio <- n_com / n_total
  ## check results
  cat(
    "GES & $\\lambda = ", lambda_use, "$&", n1, "&", n2, "&", n_com,
    "&", n_total, "&", round(n_ratio, 4), "\\\\\n"
  )
  # cat("lambda: ", lambda_use, c(sum(ges_adj1), sum(ges_adj2), sum(ges_adj)) / 2, "\n")
}
######################### joint GES method ########################
set.seed(1)
intercept_use <- FALSE
lambdas <- c(1, 2, 3, 4, 5)

## Joint GES the first step
ges_joint_fun <- function(data, lambda) {
  source("Section6/newclass.R")
  p <- ncol(data[[1]])
  dag_list <- list()
  l0score <- new("MultiGaussL0pen",
    data = data, lambda = lambda * log(p),
    intercept = intercept_use, use.cpp = FALSE
  )
  ges_fit <- ges(l0score)
  dag <- as(ges_fit$essgraph, "matrix")
  return(ifelse(dag == TRUE, 1, 0))
}

## Joint GES the second step
subset <- function(y, x, data) {
  t <- rep(0, ncol(data))
  if (length(x) <= 1) {
    t[x] <- 1
  } else {
    model <- glmnet::cv.glmnet(as.matrix(data[, x]), data[, y], family = "gaussian", intercept = intercept_use)
    nonz <- which(as.vector(coef(model)) != 0) - 1
    t[x[nonz]] <- 1
  }
  return(t)
}

# do joint estimation given single data
ges_alg <- function(dag_com, dta) {
  in_mat <- dag_com
  joint_mat <- sapply(seq_len(ncol(dta)), function(i) subset(i, which(in_mat[, i] != 0), dta))
  adj_pri <- joint_mat
  return(adj_pri)
}

for (iter in seq_len(length(lambdas))) {
  lambda_use <- lambdas[iter]
  dag_com <- ges_joint_fun(data, lambda_use)
  dag1 <- ges_alg(dag_com, data[[1]])
  dag2 <- ges_alg(dag_com, data[[2]])
  ## data set 1 results
  ges_joint_graph1 <- as(dag1, "matrix")
  ges_joint_graph1 <- ifelse(ges_joint_graph1 == 1, TRUE, FALSE)
  ges_joint_graph1 <- ges_joint_graph1 | t(ges_joint_graph1)
  n1 <- sum(ges_joint_graph1) / 2
  ## data set 2
  ges_joint_graph2 <- as(dag2, "matrix")
  ges_joint_graph2 <- ifelse(ges_joint_graph2 == 1, TRUE, FALSE)
  ges_joint_graph2 <- ges_joint_graph2 | t(ges_joint_graph2)
  n2 <- sum(ges_joint_graph2) / 2
  ## intersections
  ges_joint_graph <- ges_joint_graph1 & ges_joint_graph2
  n_com <- sum(ges_joint_graph) / 2
  n_total <- n1 + n2 - n_com
  n_ratio <- n_com / n_total
  ## check results
  cat(
    "joint GES &$\\lambda = ", lambda_use, "$&", n1, "&", n2, "&", n_com,
    "&", n_total, "&", round(n_ratio, 4), "\\\\\n"
  )
}


######################### muSuSiE method ########################
set.seed(1)
## generate graph
source("utility/graph_mcmc_multi.R")
prior_vec_list <- list()
prior_vec_list[[1]] <- c(1 / p^1.25, 1 / p^1.25, 1 / p^2)
prior_vec_list[[2]] <- c(1 / p^1.5, 1 / p^1.5, 1 / p^2.5)
prior_vec_list[[3]] <- c(1 / (2 * p^1.5), 1 / (2 * p^1.5), 1 / p^2)
prior_vec_list[[4]] <- c(1 / p^2, 1 / p^2, 1 / p^3.5)
prior_vec_list[[5]] <- c(1 / (2 * p^2), 1 / (2 * p^2), 1 / p^3.5)
prior_vec_list[[6]] <- c(1 / p^2, 1 / p^2, 1 / p^2)

scale_x <- FALSE
intercept <- TRUE # can't get two graphs if setting false
iter_max <- 1e5

#### Do MCMC with order
## get order
dta <- rbind(data[[1]], data[[2]])
# order_int <- NULL
score_ges <- new("GaussL0penObsScore", data = dta, intercept = FALSE)
ges_fit <- ges(score_ges)
ges_adj <- as(ges_fit$repr, "matrix")
ges_adj <- ifelse(ges_adj == TRUE, 1, 0)
graph_i <- igraph::graph_from_adjacency_matrix(ges_adj, mode = "directed", diag = FALSE)
order_int <- as.numeric(igraph::topo_sort(graph_i))
## Do MCMC
library(foreach)
library(doParallel)
library(doRNG)
time1 <- Sys.time()
cl <- makeCluster(length(prior_vec_list))
registerDoParallel(cl)
out_res <- foreach(iter_prior = seq_len(length(prior_vec_list))) %dorng% {
  prior_vec <- prior_vec_list[[iter_prior]]
  Graph_MCMC_multi(data,
    scale_x = scale_x, intercept = intercept,
    order_int = order_int, iter_max = iter_max, sigma02_int = NULL, sigma2_int = NULL,
    prior_vec = prior_vec, itermax = 100, tol = 1e-4, sigma0_low_bd = 1e-8,
    burn_in = iter_max - 5000
  )
}
stopCluster(cl)
print(Sys.time() - time1)
saveRDS(out_res, "Section6/results/out_mcmc.rds")
print("Finish muSuSiE-DAG method.")

## check results
for (iter_prior in seq_len(length(prior_vec_list))) {
  res_tmp <- out_res[[iter_prior]]
  prior <- prior_vec_list[[iter_prior]]
  alpha_mat_1 <- res_tmp$alpha_mat_1
  alpha_mat_2 <- res_tmp$alpha_mat_2
  A_mat_1 <- res_tmp$A_mat_1
  A_mat_2 <- res_tmp$A_mat_1
  #### Compare results
  ## data set 1
  adj_1 <- ifelse(alpha_mat_1 > 0.5, 1, 0)
  adj_1 <- t(adj_1)
  n1 <- sum(adj_1)
  ## data set 2
  adj_2 <- ifelse(alpha_mat_2 > 0.5, 1, 0)
  adj_2 <- t(adj_2)
  n2 <- sum(adj_2)
  # intersection
  n_com <- length(intersect(which(adj_1 == 1), which(adj_2 == 1)))
  n_total <- n1 + n2 - n_com
  n_ratio <- n_com / n_total
  ## check results
  cat(
    "muSuSiE-DAG &", prior[1], "&", prior[2], "&", n1, "&", n2, "&", n_com,
    "&", n_total, "&", round(n_ratio, 4), "\\\\\n"
  )
}
