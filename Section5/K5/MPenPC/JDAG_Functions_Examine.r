# functions to examine the performance of different estimation
# GraphsAccuracy: compare a sequence of estimates to the true graph
# OmegasAccuracy: compare a sequence of estimates to the true precision matrix
# CompareGraphSeqs: compare multiple sequences of estimates based on GraphsAccuracy
# AbsDiff.list: calculate absolute difference between two sequences of matrices

# Accuracy of the Graph Estimation
# FP and FN of the estimated graphs
# for a sequence of graphical models
# typically created by huge/glasso
# SHOULD be improved by considering signs of edges
GraphsAccuracy = function(graphs, graph.t, symmetrize = T) {
    if (!is.list(graphs)) {
        graphs = list(graphs)
    }
    K = length(graphs)
    p = nrow(graph.t)
    N = p * (p - 1) / ifelse(symmetrize, 2, 1)
    graph.t = as.matrix(graph.t)
    diag(graph.t) = FALSE
    if (symmetrize) { graph.t = graph.t | t(graph.t) }
    err = foreach(i = 1:K, .combine = rbind) %dopar% {
        graph.i = as.matrix(graphs[[i]])
        diag(graph.i) = FALSE
        if (symmetrize) { graph.i = graph.i | t(graph.i) }
        matrix(c(sum(graph.i & graph.t), sum(graph.i)), 1, 2)
    } / ifelse(symmetrize, 2, 1)
    err = as.data.frame(err)
    names(err) = c("TP", "DF")
    err$DF_T = sum(graph.t) / ifelse(symmetrize, 2, 1)
    err$FP = err$DF - err$TP
    err$FN = err$DF_T - err$TP
    err$FPR = err$FP / (N - err$DF_T)
    err$TPR = err$TP / err$DF_T
    try(rownames(err) <- names(graphs))
    err$Lift = with(err, N * TP/DF_T/DF)
    err[, c("TP", "FP", "FN", "DF", "DF_T", "Lift", "FPR", "TPR")]
}

# Compute the average accuracy of a bunch of graphs
MeanGraphsAccuracy = function(paths.est, paths.act, ...) {
    stopifnot(is.list(paths.act), is.list(paths.est), length(paths.act) == length(paths.est))
    K = length(paths.est)
    p = nrow(paths.act[[1]])
    # mat.i = foreach(k = 1:K, .combine = rbind) %do% {
    #     GraphsAccuracy(paths.est[[k]], paths.act[[k]], ...)
    # }
    # perf.i = colMeans(mat.i)
    # perf.i['FPR'] = sum(mat.i[,'FP']) / sum(p^2 - p - mat.i[,'DF'])
    # perf.i['TPR'] = sum(mat.i[,'TP']) / sum(mat.i[,'DF_T'])
    # return(as.data.frame(t(perf.i)))
    mat.i = foreach(k = 1:K, .combine = '+') %do% {
        as.matrix(GraphsAccuracy(paths.est[[k]], paths.act[[k]], ...))
    } / K
    mat.i = as.data.frame(mat.i)
    perf.i = mat.i
    # perf.i$FPR = mat.i$FP / (p^2 - p - mat.i$DF)
    # perf.i$TPR = mat.i$TP / mat.i$DF_T
    return(perf.i)
}

# Accuracy of the Precision Matrix Estimation
# Positive/Negative/Zero Accuracy
# Sum of Absolute Differences
# for a sequence of Omegas
OmegasAccuracy = function(Omegas.est, Omega.act){
	stopifnot(is.list(Omegas.est))
	M = length(Omegas.est)
	diag(Omega.act) = 0
	scale.act = mean(abs(Omega.act))
	
	criteria = c("AbsDiff", "SignDiff", "NormDiff", "DF", "DF_T")
	accu.mat = matrix(0, M, length(criteria))
	colnames(accu.mat) = criteria
	
	for (m in 1:M){
		Omega.m = Omegas.est[[m]]
		diag(Omega.m) = 0
		scale.m = mean(abs(Omega.m))
		
		accu.mat[m, 1] = mean(abs(Omega.m - Omega.act))
		accu.mat[m, 2] = mean(sign(Omega.m) == sign(Omega.act))
		accu.mat[m, 3] = mean(abs(Omega.m/scale.m - 
			Omega.act/scale.act))
		accu.mat[m, 4] = sum(Omega.m != 0) / 2
	}
	accu.mat[ , 5] = sum(Omega.act != 0) / 2
	return(accu.mat)
}


# depends on GraphsAccuracy
# accu.mats is a list of accuracy matrices
CompareGraphSeqs = function(accu.mats, new.window = FALSE){
	K = length(accu.mats)
	cols = rainbow(2)
	
	range.DF = numeric(2)
	range.F = numeric(2)
	for (k in 1:K){
		range.DF[1] = min(c(range.DF[1], accu.mats[[k]][ , "DF"]))
		range.DF[2] = max(c(range.DF[2], accu.mats[[k]][ , "DF"]))
		range.F[1] = min(c(range.F[1], 
			accu.mats[[k]][ , "FP"], accu.mats[[k]][ , "FN"]))
		range.F[2] = max(c(range.F[2], 
			accu.mats[[k]][ , "FP"], accu.mats[[k]][ , "FN"]))
	}
		
	if (new.window) windows()
	plot(range.DF, range.F, type = "n", xlab = "DF", ylab = "FP/FN",
		main = "False Positives/Negatives of Graph Estimation")
	for (k in 1:K){
		lines(accu.mats[[k]][ , "DF"], accu.mats[[k]][ , "FP"],
			col = cols[1], lty = k)
		lines(accu.mats[[k]][ , "DF"], accu.mats[[k]][ , "FN"],
			col = cols[2], lty = k)
	}
	legend("topleft", legend = c(paste("FP", 1:K), paste("FN", 1:K)),
		col = rep(cols, each = K), lty = rep(1:K, 2))
	
	return(0)
}

CompareTwoGraphs = function(graph.list) {
    stopifnot(length(graph.list) == 2)
    g1 = (graph.list[[1]]) | t(graph.list[[1]])
    diag(g1) = FALSE
    g2 = (graph.list[[2]]) | t(graph.list[[2]])
    diag(g2) = FALSE
    p = nrow(graph.list[[1]])
    res = data.frame(p = p,
                     DF.1 = sum(g1) / 2,
                     DF.2 = sum(g2) / 2,
                     common = sum(g1 & g2) / 2,
                     lift = mean(g1 & g2) / mean(g1) / mean(g2))
    return(res)
}

CompareMultiGraphs = function(graph.list, names = 1:length(graph.list)) {
    K = length(graph.list)
    p = nrow(graph.list[[1]])
    graph.list2 = foreach(k = 1:K) %do% {
        g.k = (graph.list[[k]]) | t(graph.list[[k]])
        diag(g.k) = FALSE
        g.k
    }
    comb.K2 = t(combn(K, 2))
    names.K2 = matrix(names[comb.K2], nrow(comb.K2), ncol(comb.K2))
    res = matrix(NA, 1, 1 + K + choose(K, 2) * 2)
    colnames(res) = c('p', paste0('df.', names), 
                      paste('common', names.K2[,1], names.K2[,2], sep = '.'),
                      paste('lift', names.K2[,1], names.K2[,2], sep = '.'))
    res[1, 'p'] = p
    for (k in 1:K) {
        res[1, paste0('df.', names[k])] = sum(graph.list2[[k]]) / 2
    }
    for (kk in 1:nrow(comb.K2)) {
        res[1, paste('common', names.K2[kk,1], names.K2[kk,2], sep = '.')] = 
            sum(graph.list2[[comb.K2[kk,1]]] & graph.list2[[comb.K2[kk,2]]]) / 2
        res[1, paste('lift', names.K2[kk,1], names.K2[kk,2], sep = '.')] = 
            mean(graph.list2[[comb.K2[kk,1]]] & graph.list2[[comb.K2[kk,2]]]) / 
            mean(graph.list2[[comb.K2[kk,1]]]) / mean(graph.list2[[comb.K2[kk,2]]])
    }
    return(as.data.frame(res))
}

mypermute = function(K, m) {
    mat = matrix(NA, choose(K, m), m)
    
}

# Comparison between two sequences of matrices
AbsDiff.list = function(A, B, inverse = FALSE){
  stopifnot(length(A) == length(B))
  K = length(A)
  if (inverse){
    for (k in 1:K){
      B[[k]] = solve(B[[k]])
    }
  }
  
  absdiff = numeric(K)
  for (k in 1:K){
    absdiff[k] = 2 * sum(abs(A[[k]] - B[[k]])) / 
      sum(abs(A[[k]]) + abs(B[[k]]))
  }
  return(absdiff)
}


