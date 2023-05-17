# Functions to Implement Hard/Soft Labeling
# RegisterLabel : register lab2 based on lab1 (canonical)
# PredMultinom  : predict based on multi-category linear classifiers
# SoftLabel.QDA : soft labeling based on hard labels and QDA
# Softlabel.NB  : soft labeling based on Naive Bayes
# SoftLabel.CIG : ... and CIG
# WSample       : generate weighted sample based on soft labeling

# Sample Expansion based on Soft Labeling

library(doParallel)
library(foreach)
library(huge)

# ExpSample = function(Y, probs) {
#     n     = nrow(Y)
#     p     = ncol(Y)
#     stopifnot(nrow(probs) == n)
#     K     = ncol(probs)
#     probs = probs / matrix(rowSums(probs), n, K, byrow = FALSE)
#     
#     Y.new = foreach(k = 1:K, .combine = rbind) %do% {
#         Y * sqrt(probs[, k])
#     }
#     return(Y.new)
# }

ExpSample = function(Y, probs, thres = 1e-8) {
    n     = nrow(Y)
    p     = ncol(Y)
    stopifnot(nrow(probs) == n)
    K     = ncol(probs)
    probs = probs / matrix(rowSums(probs), n, K, byrow = FALSE)

    Y.new = foreach(k = 1:K) %do% {
        obs.k = which(probs[,k] > thres)
        Y[obs.k,] * sqrt(probs[obs.k, k])
    }
# 	Y.new = matrix(0, nrow = n * K, ncol = p)
#     for (k in 1:K){
#       Y.new[((k-1)*n + 1):(k*n), ] = Y * sqrt(probs[,k])
#     }
    return(Y.new)
}

# weighted sample
WSample = function(Y, probs){
    n     = nrow(Y)
    p     = ncol(Y)
    stopifnot(nrow(probs) == n)
    K     = ncol(probs)
    probs = probs / matrix(rowSums(probs), n, K, byrow = FALSE)
    Y.new = list()
    S.new = list()
    
    for (k in 1:K) {
        ybar.k = c(t(probs[,k]) %*% Y) / sum(probs[,k])
        Y0.k   = Y - matrix(ybar.k, nrow = n, ncol = p, byrow = TRUE)
        S.new[[k]]= t(Y0.k) %*% diag(probs[,k]) %*% Y0.k / (sum(probs[,k]) - 1)
        Y.new[[k]] = diag(sqrt(probs[,k])) %*% Y0.k / sqrt(mean(probs[,k]))
    }
    return(list(Y = Y.new, S = S.new, weights = colMeans(probs)))
}


# ------------------------ Conditional Probs Estimation ------------------------
# For now, QDA is most convenient to estimate multi-class probs
# Output a matrix of probs
# each column corresponds to a label by label names
SoftLabel.QDA = function(dat, label) {
	label.unique = sort(unique(label))
	K            = length(label.unique)
	n            = nrow(dat)
	p            = ncol(dat)
	
	n.seq     = integer(K)
	mu.seq    = matrix(0, K, p)
	Omega.seq = list()
	
	for (k in 1:K) {
	    dat.k   = dat[label == label.unique[k], ]
		mu.k    = colMeans(dat.k)
		S.k     = var(dat.k)
		
		if (all(eigen(S.k)$values > 0)) {
		    Omega.k = solve(S.k)
		} else {
		    sgs.k   = huge(dat.k, method = "glasso", nlambda = 50, verbose = FALSE)
		    Omega.k = as.matrix(huge.select(sgs.k, criterion = "ebic")$opt.icov)
		}
	    
		n.seq[k]       = nrow(dat.k)
	    mu.seq[k, ]    = mu.k
		Omega.seq[[k]] = Omega.k
	}
	probs.mat = ProbsQDA(dat, mu.seq, Omega.seq, w = n.seq / sum(n.seq))
	if (any(is.na(probs.mat))){
		which.na = which(is.na(probs.mat), arr.ind = TRUE)
		which.na = unique(which.na[,1])
		probs.mat[which.na, ]                       = 0
		probs.mat[cbind(which.na, label[which.na])] = 1
	}
	colnames(probs.mat) = label.unique
	
	return(probs.mat)
}

# For now, QDA is most convenient to estimate multi-class probs
# Output a matrix of probs
# each column corresponds to a label by label names
SoftLabel.NB = function(dat, label) {
	label.unique = sort(unique(label))
	
	dat.xy       = data.frame(y = label, dat)
	dat.xy$y     = as.factor(dat.xy$y)
		
	fm.NB        = e1071::naiveBayes(y ~ ., data = dat.xy)
	
	probs.mat    = predict(fm.NB, dat.xy, type = "raw")
	probs.mat    = probs.mat[ ,as.character(label.unique)]
	
	colnames(probs.mat) = label.unique
	
	return(probs.mat)
}

# Predict based on a Multi-Categorical Linear Classifier
# B: (K-1)-by-(p+1) coefficient matrix; each row for a hyperplane
PredMultinom = function(B, X){
  if (ncol(matrix(X)) == 1) X = t(matrix(X))
  K = nrow(coef(fm0)) + 1
  probs.test = cbind(1, X) %*% t(coef(fm0))
  probs.test = cbind(1, exp(probs.test))
  probs.test = probs.test / matrix(rowSums(probs.test), 
                                   nrow = nrow(X), 
                                   ncol = K, byrow = FALSE)
  g.hat = (apply(probs.test, 1, which.max))
  return(list(p.hat = probs.test, g.hat = g.hat))
}


# Register lab2 based on lab1 (canonical)
# lab1 must be 1~K for prob input lab2
RegisterLabel = function(lab1, lab2, maxiter = 1e4){
	if (is.matrix(lab2)) {
	    prob.input = TRUE
	    prob2 = lab2 / rowSums(lab2)
		lab2  = apply(prob2, 1, which.max)
	} else {
	    prob.input = FALSE
	}
	unique.1 = sort(unique(lab1))
	unique.2 = sort(unique(lab2))
	stopifnot(length(unique.1) == length(unique.2), 
		length(lab1) == length(lab2))
	
	K = length(unique.1)
	n = length(lab1)
	
	stopifnot(factorial(K) <= maxiter)
	
	lab1.num = integer(n)
	lab2.num = integer(n)
	for (k in 1:K){
		lab1.num[lab1 == unique.1[k]] = k
		lab2.num[lab2 == unique.2[k]] = k
	}
	
	allperm = gtools::permutations(K, K)
	n.allperm = nrow(allperm)
	lab2.num.perm = allperm[, lab2.num]
	overlap = rowSums(lab2.num.perm ==
		matrix(lab1.num, nrow = n.allperm, ncol = n, byrow = TRUE))
	
	perm.match = which.max(overlap)[1]
	lab2.match = unique.1[lab2.num.perm[perm.match, ]]
	
	if (prob.input) {
	    prob3 = matrix(NA, nrow(prob2), ncol(prob2))
		for (k in 1:K) {
		    prob3[ , allperm[perm.match, k]] = prob2[,k]
		}
	    return(list(label = lab2.match, prob = prob3))
	} else {
	    return(lab2.match)
	}	
}

