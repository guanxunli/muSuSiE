library(pcalg)
library(parallel)
source("Section5/graph_generation.R")
p <- 100
n_tol <- 1200
K <- 5
n <- n_tol / K
e_com <- 100 # 50 100 100
e_pri <- 20 # 50 50 20

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

# ############### PC method ###########
# set.seed(2021)
# alphas <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)
# pc_fun <- function(dta, alphas = c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)) {
#   p <- ncol(dta)
#   dta_cor <- cor(dta)
#   dag_list <- list()
#   for (iter_alpha in seq_len(length(alphas))) {
#     alpha <- alphas[iter_alpha]
#     pc_fit <- pc(
#       suffStat = list(C = dta_cor, n = dim(dta)[1]),
#       indepTest = gaussCItest, alpha = alpha,
#       labels = sapply(1:p, toString)
#     )
#     dag <- as(pc_fit@graph, "matrix")
#     dag_list[[iter_alpha]] <- ifelse(dag == TRUE, 1, 0)
#   }
#   return(dag_list)
# }
# 
# cl <- makeCluster(50)
# registerDoParallel(cl)
# out_res <- foreach(iter = seq_len(n_graph)) %dorng% {
#   library(pcalg)
#   time1 <- Sys.time()
#   ## load data
#   data <- graph_sim$X[[iter]]
#   ## Do GES
#   dag_list <- list()
#   ## Do GES
#   for (iter_K in seq_len(K)) {
#     dag_list[[iter_K]] <- pc_fun(data[[iter_K]], alphas = alphas)
#   }
#   time_use <- Sys.time() - time1
#   return(list(dag_list = dag_list, time_use = time_use))
# }
# stopCluster(cl)
# saveRDS(out_res, paste0("Section5/K5/results/pc_com", e_com, "pri", e_pri, ".rds"))
# print("Finish PC.")
# 
# ############### GES method ###########
# set.seed(2021)
# lambdas <- c(1, 2, 3, 4, 5)
# ges_fun <- function(dta, lambdas = c(1, 2, 3, 4, 5)) {
#   p <- ncol(dta)
#   dag_list <- list()
#   for (iter_lambda in seq_len(length(lambdas))) {
#     lambda <- lambdas[iter_lambda]
#     l0score <- new("GaussL0penObsScore", data = dta, lambda = lambda * log(p), intercept = FALSE)
#     ges_fit <- ges(l0score)
#     dag <- as(ges_fit$repr, "matrix")
#     dag_list[[iter_lambda]] <- ifelse(dag == TRUE, 1, 0)
#   }
#   return(dag_list)
# }
# 
# cl <- makeCluster(50)
# registerDoParallel(cl)
# out_res <- foreach(iter = seq_len(n_graph)) %dorng% {
#   time1 <- Sys.time()
#   library(pcalg)
#   ## load data
#   data <- graph_sim$X[[iter]]
#   dag_list <- list()
#   ## Do GES
#   for (iter_K in seq_len(K)) {
#     dag_list[[iter_K]] <- ges_fun(data[[iter_K]], lambdas = lambdas)
#   }
#   time_use <- Sys.time() - time1
#   return(list(dag_list = dag_list, time_use = time_use))
# }
# stopCluster(cl)
# saveRDS(out_res, paste0("Section5/K5/results/GES_com", e_com, "pri", e_pri, ".rds"))
# print("Finish GSE.")
# 
# ############### joint GES ###########
# set.seed(2021)
# #### joint GES method the first step
# lambdas <- c(1, 2, 3, 4, 5)
# ges_joint_fun <- function(data, lambdas = c(1, 2, 3, 4, 5)) {
#   source("Section5/newclass.R")
#   p <- ncol(data[[1]])
#   dag_list <- list()
#   for (iter_lambda in seq_len(length(lambdas))) {
#     lambda <- lambdas[iter_lambda]
#     l0score <- new("MultiGaussL0pen",
#       data = data, lambda = lambda * log(p),
#       intercept = FALSE, use.cpp = FALSE
#     )
#     ges_fit <- ges(l0score)
#     dag <- as(ges_fit$essgraph, "matrix")
#     dag_list[[iter_lambda]] <- ifelse(dag == TRUE, 1, 0)
#   }
#   return(dag_list)
# }
# 
# ## Joint GES the second step
# subset <- function(y, x, data) {
#   t <- rep(0, ncol(data))
#   if (length(x) <= 1) {
#     t[x] <- 1
#   } else {
#     model <- glmnet::cv.glmnet(as.matrix(data[, x]), data[, y], family = "gaussian", intercept = FALSE)
#     nonz <- which(as.vector(coef(model)) != 0) - 1
#     t[x[nonz]] <- 1
#   }
#   return(t)
# }
# 
# # do joint estimation given single data
# ges_alg <- function(dag_list, dta) {
#   adj_list <- list()
#   for (iter in seq_len(length(dag_list))) {
#     in_mat <- dag_list[[iter]]
#     joint_mat <- sapply(seq_len(ncol(dta)), function(i) subset(i, which(in_mat[, i] != 0), dta))
#     adj_list[[iter]] <- joint_mat
#   }
#   return(adj_list)
# }
# 
# cl <- makeCluster(50)
# registerDoParallel(cl)
# out_res <- foreach(iter = seq_len(n_graph)) %dorng% {
#   library(pcalg)
#   time1 <- Sys.time()
#   ## load data
#   data <- graph_sim$X[[iter]]
#   dag_list <- list()
#   ## Do joint GES the first step
#   dag_list_com <- ges_joint_fun(data)
#   ## Do joint GES the second step
#   for (iter_K in seq_len(K)) {
#     dag_list[[iter_K]] <- ges_alg(dag_list_com, data[[iter_K]])
#   }
#   time_use <- Sys.time() - time1
#   return(list(dag_list = dag_list, time_use = time_use))
# }
# stopCluster(cl)
# saveRDS(out_res, paste0("Section5/K5/results/jointGES_com", e_com, "pri", e_pri, ".rds"))
# print("Finish joint GSE.")
# 
# ############### muSuSiE ###########
# source("utility/graph_mcmc_multi.R")
# #### define prior
# prior_pi_list <- list()
# prior_pi_list[[1]] <- c(1 / (5 * p^1.5), 1 / (10 * p^1.75), 1 / (10 * p^2), 1 / (5 * p^2.25), 1 / p^2.5)
# prior_pi_list[[2]] <- c(1 / p^1.75, 1 / p^2, 1 / p^2.25, 1 / p^2.5, 1 / p^3)
# prior_pi_list[[3]] <- c(1 / (5 * p^1.75), 1 / (10 * p^2), 1 / (10 * p^2.25), 1 / (5 * p^2.5), 1 / p^3)
# prior_pi_list[[4]] <- c(1 / p^2, 1 / p^2.25, 1 / p^2.5, 1 / p^2.75, 1 / p^3)
# n_prior <- length(prior_pi_list)
# prior_vec_list <- list()
# for (iter_pi in seq_len(n_prior)) {
#   n_group <- 2^K - 1
#   com_list <- list()
#   com_mat <- matrix(c(0, 1), ncol = 1)
#   for (iter in 2:K) {
#     com_mat_copy <- com_mat
#     com_mat <- cbind(1, com_mat)
#     com_mat_copy <- cbind(0, com_mat_copy)
#     com_mat <- rbind(com_mat_copy, com_mat)
#   }
#   com_mat <- com_mat[-1, ]
# 
#   for (iter_com in seq_len(n_group)) {
#     com_list[[iter_com]] <- which(com_mat[iter_com, ] == 1)
#   }
#   com_length <- lapply(com_list, length)
#   prior_vec <- rep(NA, n_group)
#   prior_pi <- prior_pi_list[[iter_pi]]
#   for (iter in seq_len(n_group)) {
#     prior_vec[iter] <- prior_pi[com_length[[iter]]]
#   }
#   prior_vec_list[[iter_pi]] <- prior_vec
# }
# #### simulations
# iter_max <- 1e5
# for (iter_prior in seq_len(n_prior)) {
#   set.seed(2021)
#   print(iter_prior)
#   out_res <- list()
#   prior_vec <- prior_vec_list[[iter_prior]]
#   ## do parallel
#   cl <- makeCluster(50)
#   registerDoParallel(cl)
#   out_res <- foreach(iter = seq_len(n_graph)) %dorng% {
#     library(pcalg)
#     time1 <- Sys.time()
#     # get order
#     dta <- matrix(NA, nrow = K * n, ncol = p)
#     adj_true <- list()
#     for (iter_K in seq_len(K)) {
#       dta[(1 + (iter_K - 1) * n):(iter_K * n), ] <- graph_sim$X[[iter]][[iter_K]]
#       adj_true[[iter_K]] <- t(graph_sim$G[[iter]][[iter_K]])
#     }
#     score_ges <- new("GaussL0penObsScore", data = dta, intercept = FALSE)
#     ges_fit <- ges(score_ges)
#     ges_adj <- as(ges_fit$repr, "matrix")
#     ges_adj <- ifelse(ges_adj == TRUE, 1, 0)
#     graph_i <- igraph::graph_from_adjacency_matrix(ges_adj, mode = "directed", diag = FALSE)
#     order_int <- as.numeric(igraph::topo_sort(graph_i))
#     # Do MCMC
#     res <- Graph_MCMC_multi(
#       dta_list = graph_sim$X[[iter]], scale_x = FALSE, intercept = TRUE,
#       order_int = order_int, iter_max = iter_max,
#       prior_vec = prior_vec, itermax = 100, L_max = 10,
#       burn_in = iter_max - 5000
#     )
#     time_use <- Sys.time() - time1
#     return(list(res = res, time_use = time_use))
#   }
#   stopCluster(cl)
#   saveRDS(out_res, paste0("Section5/K5/results/muSuSiE_prior", iter_prior, "com", e_com, "pri", e_pri, ".rds"))
# }

############### JESC method ###########
set.seed(2021)
cl <- makeCluster(50)
registerDoParallel(cl)
out_res <- foreach(iter = seq_len(n_graph)) %dorng% {
  source("Section5/K5/utility_esc.R")
  # hyperparameters
  alpha = 0.999
  gamma = 1
  nu0 = 0.1
  c1 = 1 #c1 in ESC paper, doesn't matter when set to 1
  c2 = 2
  b = 1/(p*(K - 1)) #c2 in the paper
  niter = 10000
  nburn = 0.2 * niter
  nadap = 0
  time1 <- Sys.time()
  ## load data
  data <- graph_sim$X[[iter]]
  X = array(0, dim = c(K, dim(data[[1]])))
  for (iter_K in seq_len(K)) {
    X[iter_K, , ] <- data[[iter_K]]
  }
  ## run JESC
  res = mdag(X, alpha, gamma, nu0, c1, c2, c3 = NULL, b, niter, nburn)
  ## sort return list
  l1 = list()
  for (k in 1:K) {
    l2 = list(NULL)
    for (j in 1:(p-1)) {
      m <- matrix(0, niter, j)
      for (t in (nburn + 1):(nburn + niter)) {
        m[t-nburn, ] = res[[j]][[t]][k, ]
      }
      l2[[j]] = m
    }
    l1[[k]] = l2
  }
  ## generate DAG
  dag_list <- list()
  for (iter_K in seq_len(K)) {
    dag_tmp <- matrix(0,p,p)
    for(j in 2:p){
      Sj_mat = l1[[iter_K]][[j-1]]
      dag_tmp[j,seq(1:(j-1))] <- apply(Sj_mat, 2, mean)
    }
    dag_list[[iter_K]] <- dag_tmp
  }
  time_use <- Sys.time() - time1
  return(list(dag_list = dag_list, time_use = time_use))
}
stopCluster(cl)
saveRDS(out_res, paste0("Section5/K5/results/JESC_com", e_com, "pri", e_pri, ".rds"))
print("Finish JESC method.")