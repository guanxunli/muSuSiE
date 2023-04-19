# K : number of graphs
# n_tol : total number of observations
# p : number of nodes
# e_com : number of common edges
# e_pri : number of private edges
# n_graph : number of simulations
# min_sig : minimun signal of the weights
# low_err_var : lower bound for error variance
# upp_err_var : upper bound for the error variance
graph_generation <- function(K = 2, n_tol = 600, p = 100, e_com = 100, e_pri = 30, n_graph = 1, min_sig = 0.1,
                             low_err_var = 1, upp_err_var = 2.25) {
  ## initialization
  n <- n_tol / K
  # G is the true graph (p x p)
  # X is the data set from the true graph (n x p)
  # A is the coefficients matrix (p x p)
  X <- list()
  G <- list()
  A <- list()
  ## begin iteration
  for (iter_graph in seq_len(n_graph)) {
    G[[iter_graph]] <- list()
    A[[iter_graph]] <- list()
    X[[iter_graph]] <- list()
    edge_com <- sample(seq_len(p * (p - 1) / 2), e_com)
    for (iter_K in seq_len(K)) {
      # generate graph
      G[[iter_graph]][[iter_K]] <- matrix(0, nrow = p, ncol = p)
      edge_pri <- sample(setdiff(seq_len(p * (p - 1) / 2), edge_com), e_pri)
      G[[iter_graph]][[iter_K]][lower.tri(G[[iter_graph]][[iter_K]])][c(edge_com, edge_pri)] <- 1
      # generate weight
      A[[iter_graph]][[iter_K]] <- matrix(0, nrow = p, ncol = p)
      A[[iter_graph]][[iter_K]][which(G[[iter_graph]][[iter_K]] == 1)] <- 
        2 *  (rbinom(e_com + e_pri, 1, 0.5) - 1/2) * runif(e_com + e_pri, min = min_sig, max = 1)
      # generate data
      err_var <- runif(1, min = low_err_var, max = upp_err_var)
      err_vec <- matrix(rnorm(n * p, mean = 0, sd = sqrt(err_var)), ncol = n)
      X[[iter_graph]][[iter_K]] <- solve(diag(1, nrow = p) - A[[iter_graph]][[iter_K]], err_vec)
      X[[iter_graph]][[iter_K]] <- t(X[[iter_graph]][[iter_K]])
    }
  }
  # return graph
  return(list(G = G, A = A, X = X))
}