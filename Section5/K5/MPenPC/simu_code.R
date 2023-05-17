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

#### MPenPC method ####
setwd("Section5/K5/MPenPC/")
source("JDAG_head.r")
cl <- makeCluster(50)
registerDoParallel(cl)
out_res <- foreach(iter = seq_len(n_graph)) %dorng% {
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
  dag_list <- res$graphs1
  time_use <- Sys.time() - time1
  return(list(dag_list = dag_list, time_use = time_use))
}
stopCluster(cl)
saveRDS(out_res, paste0("Section5/K5/results/MPenPC_com", e_com, "pri", e_pri, ".rds"))
print("Finish MPenPC")