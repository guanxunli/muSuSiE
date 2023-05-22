library(pcalg)
library(parallel)
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

#### MPenPC method ####
setwd("Section5/K2/MPenPC/")
cl <- makeCluster(50)
registerDoParallel(cl)
out_res <- foreach(iter = seq_len(n_graph)) %dorng% {
  source("JDAG_head.r")
  time1 <- Sys.time()
  ## load data
  data <- graph_sim$X[[iter]]
  dta_use <- matrix(NA, nrow = n_tol, ncol = p)
  for (iter_K in seq_len(K)) {
    dta_use[(1 + n * (iter_K - 1)):(n * iter_K), ] <- data[[iter_K]]
  }
  probs_mat <- matrix(0, nrow = n_tol, ncol = K)
  for (iter_K in seq_len(K)) {
    probs_mat[(1 + n * (iter_K - 1)):(n * iter_K), iter_K] <- 1
  }
  res <- SoftGraphs(dat = dta_use, probs.mat = probs_mat, joint = T)
  res_list <- list(res)
  res_list <- lapply(res_list, SymmetrizeRes, 'Union')
  res_list <- lapply(res_list, function(x) list(graphs = x$graphs1, probs = x$probs))
  dag_list <- SoftPC(res_list[[1]], X = dta_use, alpha = 0.02)
  dag_list <- lapply(dag_list, function(x) as(x@graph, 'matrix'))
  time_use <- Sys.time() - time1
  return(list(dag_list = dag_list, time_use = time_use))
}
stopCluster(cl)
saveRDS(out_res, paste0("../results/MPenPC_com", e_com, "pri", e_pri, ".rds"))
print("Finish MPenPC")
