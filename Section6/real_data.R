# # load data
load("Section6/ovarian.rda")
p <- ncol(data[[1]])
library(pcalg)
library(graph)
library(parallel)
library(stabs)

################################ PC method ########################
set.seed(1)
cutoff_vec <- seq(0.6, 0.9, by = 0.05)
## PC input
stab_input <- function(i) {
  p <- ncol(data[[i]])
  dt <- cbind(as.matrix(data[[i]]), matrix(0, nrow = nrow(data[[i]]), ncol = p * (p - 1)))
  return(list(x = dt, y = rep(i, nrow(data[[i]]))))
}
stab_input_list <- lapply(seq_len(length(data)), stab_input)

## learn causal networks
stabs_pc <- function(x, y, q, ...) {
  # Y is the label of the classes, X is the input matrix
  idx <- y[1]
  totcol <- nrow(data[[idx]])
  dt <- data[[idx]][sample(1:totcol, as.integer(0.9 * totcol), replace = FALSE), ]
  # dt <- data[[idx]]
  p <- ncol(dt)

  # train the model
  alphas <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)
  model_alpha <- function(alpha) {
    pc_fit <- pc(
      suffStat = list(C = cor(dt), n = nrow(dt)),
      indepTest = gaussCItest, alpha = alpha,
      labels = sapply(1:p, toString)
    )
    dag <- as(pc_fit@graph, "matrix")
    as.vector(dag != 0)
  }

  # get the path and selected variables
  path <- sapply(alphas, model_alpha)
  selected <- rowSums(path) != 0
  return(list(selected = selected, path = path))
}

pcdag_fun <- function(cutoff) {
  res <- lapply(
    stab_input_list,
    function(stab_input) stabsel(x = stab_input$x, y = stab_input$y, fitfun = stabs_pc, cutoff = cutoff, PFER = 1)
  )
  return(res)
}

time1 <- Sys.time()
pcdag_list <- mclapply(cutoff_vec, pcdag_fun,
  mc.cores = length(cutoff_vec)
)
print(Sys.time() - time1)
saveRDS(pcdag_list, "Section6/results/out_pc.rds")
print("Finish PC method.")

cutoff_vec2 <- c(0.55, 0.75)
for (iter in seq_len(length(cutoff_vec))) {
  pcdag_list_tmp <- pcdag_list[[iter]]
  for (iter2 in seq_len(length(cutoff_vec2))) {
    cutoff <- cutoff_vec2[iter2]
    ## data set 1 results
    pc_adj1 <- matrix(as.vector(pcdag_list_tmp[[1]]$max > cutoff), nrow = p, ncol = p)
    pc_adj1 <- pc_adj1 | t(pc_adj1)
    n1 <- sum(pc_adj1) / 2
    ## data set 2
    pc_adj2 <- matrix(as.vector(pcdag_list_tmp[[2]]$max > cutoff), nrow = p, ncol = p)
    pc_adj2 <- pc_adj2 | t(pc_adj2)
    n2 <- sum(pc_adj2) / 2
    ## intersections
    pc_adj <- pc_adj1 & pc_adj2
    n_com <- sum(pc_adj) / 2
    n_total <- n1 + n2 - n_com
    n_ratio <- n_com / n_total
    cat(
      "PC &", cutoff_vec[iter], "&", cutoff, "&", n1, "&", n2, "&", n_com,
      "&", n_total, "&", round(n_ratio, 4), "\\\\\n"
    )
  }
}

######################### GES method ########################
set.seed(1)
intercept_use <- FALSE
cutoff_vec <- seq(0.6, 0.9, by = 0.05)

## GES input
stab_input <- function(i) {
  p <- ncol(data[[i]])
  dt <- cbind(as.matrix(data[[i]]), matrix(0, nrow = nrow(data[[i]]), ncol = p * (p - 1)))
  return(list(x = dt, y = rep(i, nrow(data[[i]]))))
}
stab_input_list <- lapply(seq_len(length(data)), stab_input)

# stable GES function
stab_ges <- function(x, y, q, ...) {
  # Y is the label of the classes, X is the input matrix
  idx <- y[1]
  totcol <- nrow(data[[idx]])
  dt <- data[[idx]][sample(1:totcol, as.integer(0.9 * totcol), replace = FALSE), ]
  # dt <- data[[idx]]

  # train the model
  lambdas <- c(1, 2, 3, 4, 5)
  model_lambda <- function(lambda) {
    l0score <- new("GaussL0penObsScore", data = dt, lambda = lambda * log(ncol(dt)), intercept = intercept_use)
    ges_fit <- ges(l0score)
    dag <- as(ges_fit$essgraph, "matrix")
    as.vector(dag != 0)
  }

  # get the path and selected variables
  path <- sapply(lambdas, model_lambda)
  selected <- rowSums(path) != 0
  return(list(selected = selected, path = path))
}

gesdag_fun <- function(cutoff) {
  return(lapply(
    stab_input_list,
    function(stab_input) stabsel(x = stab_input$x, y = stab_input$y, fitfun = stab_ges, cutoff = cutoff, PFER = 1)
  ))
}
time1 <- Sys.time()
gesdag_list <- mclapply(cutoff_vec, gesdag_fun,
  mc.cores = length(cutoff_vec)
)
saveRDS(gesdag_list, "Section6/results/out_ges.rds")
print(Sys.time() - time1)
print("Finish GES method.")

#### check results
cutoff_vec2 <- c(0.55, 0.75)
for (iter in seq_len(length(cutoff_vec))) {
  gesdag_list_tmp <- gesdag_list[[iter]]
  for (iter2 in seq_len(length(cutoff_vec2))) {
    cutoff <- cutoff_vec2[iter2]
    ## data set 1 results
    ges_adj1 <- matrix(as.vector(gesdag_list_tmp[[1]]$max > cutoff), nrow = p, ncol = p)
    ges_adj1 <- ges_adj1 | t(ges_adj1)
    n1 <- sum(ges_adj1) / 2
    ## data set 2
    ges_adj2 <- matrix(as.vector(gesdag_list_tmp[[2]]$max > cutoff), nrow = p, ncol = p)
    ges_adj2 <- ges_adj2 | t(ges_adj2)
    n2 <- sum(ges_adj2) / 2
    ## intersections
    ges_adj <- ges_adj1 & ges_adj2
    n_com <- sum(ges_adj) / 2
    n_total <- n1 + n2 - n_com
    n_ratio <- n_com / n_total
    cat(
      "GES &", cutoff_vec[iter], "&", cutoff, "&", n1, "&", n2, "&", n_com,
      "&", n_total, "&", round(n_ratio, 4), "\\\\\n"
    )
  }
}

######################### joint GES method ########################
set.seed(1)
intercept_use <- FALSE
source("Section6/newclass.R")
cutoff_vec <- seq(0.6, 0.9, by = 0.05)

## learn causal networks
stabs_ges <- function(x, y, q, ...) {
  sample_data <- function(sing_dt) {
    totcol <- nrow(sing_dt)
    sing_dt[sample(1:totcol, as.integer(0.9 * totcol), replace = FALSE), ]
  }
  # Y is the label of the classes, X is the input matrix
  dt <- lapply(data, sample_data)
  # dt <- data
  lambdas <- c(1, 2, 3, 4, 5)
  model_lambda <- function(lambda) {
    l0score <- new("MultiGaussL0pen", data = dt, lambda = lambda * log(ncol(dt[[1]])), intercept = intercept_use, use.cpp = FALSE)
    ges_fit <- ges(l0score)
    dag <- as(ges_fit$essgraph, "matrix")
    as.vector(dag != 0)
  }
  path <- sapply(lambdas, model_lambda)
  selected <- rowSums(path) != 0
  return(list(selected = selected, path = path))
}

## joint GES the first step
# construct x
x <- do.call(rbind, data)
x <- cbind(x, matrix(0, nrow = nrow(x), ncol = p * (p - 1)))
# construct y
y <- c()
for (i in seq_len(length(data))) {
  y <- c(y, rep(i, nrow(data[[i]])))
}

## stable joint GES
stab_fun <- function(cutoff) {
  return(stabsel(x = x, y = y, fitfun = stabs_ges, cutoff = cutoff, PFER = 1))
}
time1 <- Sys.time()
stab_result_list <- mclapply(cutoff_vec, stab_fun, mc.cores = length(cutoff_vec))
print(Sys.time() - time1)
saveRDS(stab_result_list, "Section6/results/out_ges_joint.rds")
print("Finish the 1st step of the joint GES.")

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
ges_alg <- function(data, dag) {
  in_mat <- as(pdag2dag(dag)$graph, "matrix")
  joint_mat <- lapply(data, function(dt) sapply(seq_len(ncol(dt)), function(i) subset(i, which(in_mat[, i] != 0), dt)))
  return(lapply(joint_mat, function(sing_mat) dag2cpdag(as(sing_mat, "graphNEL"))))
}

cutoff_vec2 <- c(0.55, 0.75)
## check results
time1 <- Sys.time()
for (iter in seq_len(length(cutoff_vec))) {
  stab_result <- stab_result_list[[iter]]
  for (iter2 in seq_len(length(cutoff_vec2))) {
    cutoff <- cutoff_vec2[iter2]
    dag <- matrix(as.vector(stab_result$max > cutoff), nrow = p, ncol = p)
    dag <- as(dag, "graphNEL")

    gesdag <- ges_alg(data, dag)
    ## data set 1 results
    ges_joint_graph1 <- gesdag[[1]]
    ges_joint_graph1 <- as(ges_joint_graph1, "matrix")
    ges_joint_graph1 <- ifelse(ges_joint_graph1 == 1, TRUE, FALSE)
    ges_joint_graph1 <- ges_joint_graph1 | t(ges_joint_graph1)
    n1 <- sum(ges_joint_graph1) / 2

    ## data set 2
    ges_joint_graph2 <- gesdag[[2]]
    ges_joint_graph2 <- as(ges_joint_graph2, "matrix")
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
      "joint GES &", cutoff_vec[iter], "&", cutoff, "&", n1, "&", n2, "&", n_com,
      "&", n_total, "&", round(n_ratio, 4), "\\\\\n"
    )
  }
}
print(Sys.time() - time1)