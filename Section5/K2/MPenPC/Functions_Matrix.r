# Frequently used functions for matrix manipulation
# logdet: log of determinant
# cov2cor: transform a covariance matrix to a correlation one


# return the log determinant of a square matrix X
logdet = function(X){
	ld = determinant(X, logarithm = TRUE)
	if (ld$sign == 1){
		return(as.numeric(ld$modulus))
	}else{
		warning("Negative determinant! NA returned.")
		return(NA)
	}
}


# Transform a covariance matrix to a correlation matrix
cov2cor = function(Sigma, tol = 1e-5){
	stopifnot(isSymmetric(Sigma, tol = mean(abs(Sigma)) * tol))
	# eig.S = eigen(Sigma)
	# stopifnot(all(eig.S$values >= 0))
	
	p = nrow(Sigma)
	d = diag(Sigma)
	ix.nonc = which(d > 0)
	D = diag(1 / sqrt(d[ix.nonc]))
	
	cor.mat = matrix(0, p, p)
	cor.mat[ix.nonc, ix.nonc] = D %*% Sigma[ix.nonc, ix.nonc] %*% D
	
	return(cor.mat)
}

