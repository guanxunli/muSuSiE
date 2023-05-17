# Functions for Gaussian Sampling and Calculation
# SimGaussian: generate sample from mixture Gaussian
# ProbsQDA:    calculate class probabilities of Gaussian

library(mvtnorm)


# Generate sample from a mixture Gaussian based on PRECISION matrices
SimGaussian.O = function(n.vec, mus, Omegas) {
	stopifnot(length(n.vec) == nrow(mus), nrow(mus) == length(Omegas))
	p = ncol(mus)
	K = nrow(mus)
	n = sum(n.vec)
	X = matrix(0, nrow = n, ncol = p)
	y = integer(n)
	S = list()
	for (k in 1:K){
	    eigen.k = eigen(Omegas[[k]])
		Sigma.k = with(eigen.k, vectors %*% diag(1 / values) %*% t(vectors))
		# Sigma.k = AsSymmetric(solve(Omegas[[k]]), method = "average")
		eig.k = eigen(Sigma.k)
		if (any(eig.k$values < 0)){
			warning(paste0("Minumum eigenvalue of Omegas[[", k, "]] is ", min(eig.k$values)))
			Sigma.k = Sigma.k - min(eig.k$values) * diag(p)
		}
		ix.k = (sum(n.vec[0:(k-1)]) + 1):sum(n.vec[1:k])
		y[ix.k] = k
		X[ix.k, ] = rmvnorm(n.vec[k], mus[k, ], Sigma.k)
		S[[k]] = cov(X[ix.k, ])
	}
	rand.ix = sample(n, n, replace = FALSE)
	return(list(X = X[rand.ix, ], g = y[rand.ix], S = S))
}

# Generate sample from a mixture Gaussian based on structure equations
SimGaussian.B = function(n.vec, nus, Bs, err = 1) {
    stopifnot(length(n.vec) == nrow(nus), nrow(nus) == length(Bs))
    p = ncol(nus)
    K = nrow(nus)
    n = sum(n.vec)
    X = matrix(0, nrow = n, ncol = p)
    y = integer(n)
    S = list()
    for (k in 1:K){
        B.k = -Bs[[k]]
        B.k[upper.tri(B.k, diag = TRUE)] = 0
        ix.k = (sum(n.vec[0:(k-1)]) + 1):sum(n.vec[1:k])
        y[ix.k] = k
        X[ix.k, 1] = nus[k, 1] + rnorm(n.vec[k]) * err
        for (j in 2:p) {
            X[ix.k, j] = nus[k, j] + X[ix.k, 1:j] %*% B.k[j, 1:j] + rnorm(n.vec[k]) * err
        }
        S[[k]] = cov(X[ix.k, ])
    }
    rand.ix = sample(n, n, replace = FALSE)
    return(list(X = X[rand.ix, ], g = y[rand.ix], S = S))
}

# X = train$X[,ix.uni]; mus = mu.all[,ix.uni]; Omegas = lapply(Omegas, function(x) x[ix.uni,ix.uni]); w = NULL
# Calculate class probabilities based on given parameter settings
ProbsQDA = function(X, mus, Omegas, w = NULL){
	stopifnot(nrow(mus) == length(Omegas))
	K = nrow(mus)
	n = nrow(X)
	probs = matrix(1, nrow = n, ncol = K)
	odds = matrix(0, nrow = n, ncol = K)
	
	if (is.null(w)){
		w = rep(1, K)
	}
	if (K > 1) {
	    for (k in 2:K){
	        odds[,k] = log(w[k]) - log(w[1]) + 
	            ldmvnorm(X, mus[k, ],Omegas[[k]]) - 
	            ldmvnorm(X, mus[1, ], Omegas[[1]])
	    }
    }
	probs = exp(odds)
	if (any(is.na(probs) | (probs == Inf))){
		which.na = which(is.na(probs) | (probs == Inf), arr.ind = TRUE)
		which.na = unique(which.na[ ,1])
		label.na = apply(odds[which.na, ], 1, which.max)
		probs[which.na, ] = 0
		probs[cbind(which.na, label.na)] = 1
	}
	probs = probs / rowSums(probs)
	return(probs)
}

ldmvnorm = function(X, mu, Omega){
	n = nrow(X)
	p = ncol(X)
	mu.mat = matrix(mu, nrow = n, ncol = p, byrow = TRUE)
	d = -p/2 * log(2*pi) + logdet(Omega) / 2 -
		diag((X - mu.mat) %*% Omega %*% t(X - mu.mat)) / 2
	return(d)
}
