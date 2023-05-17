# Functions to implement the PC step of multiPenPC

# library(ParallelPC)
SoftPC = function(res.ini, X, alpha = 0.01, parallel = FALSE, n.workers = 1, PenPC = FALSE, ...) {
    n      = nrow(X)
    p      = ncol(X)
    K      = length(res.ini$graphs)
	res.ini$probs = res.ini$probs / rowSums(res.ini$probs)
    
    cors   = foreach(k = 1:K) %do% {
                 cor(diag(res.ini$probs[,k]) %*% X)
             }
    
    res.PC = list()
    for (k in 1:K) {
        suffStat.k   = list(C = cors[[k]], n = n)
        gaps.k       = as.matrix(res.ini$graphs[[k]] == 0)
        attr(gaps.k, 'dimnames') = NULL
        diag(gaps.k) = TRUE
        edges.k      = (!gaps.k) + 0
        
        
        if (!PenPC) {
            if (parallel) {
                pc.k = pc_parallel(suffStat.k, fixedGaps = gaps.k,
                                   indepTest = pcalg::gaussCItest, p = as.integer(p), 
                                   skel.method = "parallel", alpha = alpha, num_workers = n.workers, ...)   
            } else {
                pc.k = pc_stable(suffStat.k, fixedGaps = gaps.k,
                                 indepTest = pcalg::gaussCItest, p = as.integer(p), 
                                 skel.method = "stable", alpha = alpha, ...)
            }
        } else {
            pc.k = skeletonPENstable(suffStat.k, edgeWeights = edges.k,
                             indepTest = pcalg::gaussCItest, p = as.integer(p), 
                             alpha = alpha, ...)
        }
        res.PC[[k]]  = pc.k
    }
    
    return(res.PC)
}


SoftCor = function(probs, X) {
    n = nrow(X)
    K = ncol(probs)
    stopifnot(nrow(probs) == n)
    
    cors = foreach(k = 1:K) %do% {
        X.k = diag(probs[,k]) %*% X
        cor(X.k)
    }
    
    return(cors)
}

graph2Matrix = function(graph, symmetric = TRUE) {
    if (class(graph)[1] == 'graphNEL') {
        graph = igraph::graph_from_graphnel(graph)
    } else if (class(graph)[1] == 'igraph') {
        graph = igraph::graph_from_graphnel(igraph::as_graphnel(graph))
    } else {
        stop('Unsupported graph class: ', class(graph))
    }
    Mat = igraph::as_adjacency_matrix(graph)
    return(Mat)
}
