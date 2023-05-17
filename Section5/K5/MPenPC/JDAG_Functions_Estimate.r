# Functions for Joint DAG Estimation
# SoftGraphs: estimate multiple graphs based on soft labeled data
# WJGL      : weighted JGL based on soft labeled data

# Note: there are some issues with the functions, mainly on the intercept term of grpreg.
# grpreg automatically add an intercept term to the design matrix, which is not desired
# in our case since each observation has a specific weight. The intercept term shall be
# created proportional to the square root of the weight vector, regardless of separate or
# common intercept model. 
# Similar issues exist in EBIC calculation.

library(Matrix)

# ------------------ Graphical Model Estimation -------------------

# -------------------------- Test Area ----------------------------
# dat = train$X; label = NULL; probs.mat = slab.NB; joint = TRUE; penalty = "gel"; n.lambda = n.graph; n.tau = 5; intercept.common = TRUE
# dat = train$X; label = NULL; penalty = "gel"; nlambda = n.graph; ntau = 11
# -------------------------- Test Area ----------------------------
# Both XX and yy have been extended and weighted
# intercept.common: XX contains K intercept columns if FALSE
# var.common: Treat them as K * n.lambda models if TRUE
EBIC.wX = function(XX, yy, B, probs, gamma.ebic = NULL, group = NULL, 
                   var.common = TRUE, intercept.common = TRUE) {
    n.lambda = ncol(B)
    K        = ncol(probs)
    n        = nrow(probs)
    p        = ncol(XX) / K - 1 + intercept.common
    ms       = colSums(probs)
    stopifnot(nrow(B) == 1 + ncol(XX))
    if (is.null(gamma.ebic)) {
        if (var.common) {
            gamma.ebic = 1 - log(n) / log(p*K) / 2
        } else {
            gamma.ebic = 1 - log(ms) / log(p) / 2
        }
    }
    YY       = matrix(yy, length(yy), n.lambda, byrow = FALSE)
    Err2     = (YY - cbind(1, XX) %*% B)^2
    B1       = B[-1, ]
    if (intercept.common) {
        B1 = B1[-seq(1, K*(p+1), by = p+1), ]
    }
    nz.B = foreach(k = 1:K, .combine = c) %do% {
        colSums(B[(p*(k-1)+1):(p*k), ] != 0)
    }
    nz.B = matrix(nz.B, K, n.lambda)
    if (var.common) {
        nz.B = colSums(nz.B)
        s2   = colSums(Err2) / n
        res  = (log(2*pi)+1) * n + n * log(s2) + nz.B * log(n) + 
            2 * gamma.ebic * lchoose(p*K, nz.B)
    } else {
        s2  = foreach(k = 1:K, .combine = rbind) %do% {
            colSums(Err2[group[[k]], ]) / ms[k] # error
        }
        res = (log(2*pi)+1) * sum(ms) + t(ms) %*% log(s2) + t(log(ms)) %*% nz.B + 
            2 * t(gamma.ebic) %*% lchoose(p, nz.B)
    }
    
    return(c(res))
}

# Joint Estimation with Soft Label
# var.common: whether the groups share the error variance, only affect the EBIC computation
# intercept.common: whether the groups share the intercept, affect the model fitting
# SoftGraphs = function(dat, label = NULL, probs.mat = NULL, joint = F, penalty = "gel",
#                       n.lambda = 20, lambda.min = 0.05, verbose = F) {
#     # Initialization
#     K       = ncol(probs.mat)
#     dat     = as.matrix(dat)
#     n       = nrow(dat)
#     p       = ncol(dat)
#     taus    = seq(0.1, 1, by = 0.1)
#     n.tau   = length(taus)
#     probs.mat = probs.mat / matrix(rowSums(probs.mat), n, K, byrow = FALSE)
#     dat.w   = scale(ExpSample(dat, probs.mat), scale = F)
#     dat.w   = foreach(k = 1:K) %do% dat.w[((k-1)*n+1):(k*n), ]
#     
#     if (joint) {
#         group.label = rep(1:(p-1), K)
#     } else {
#         group.label = 1:(p-1)
#     }
# 
#     # Nodewise regression (jointly)
#     EBIC.wX = EBIC.wX
#     est.all = foreach(j = 1:p, .packages = c('grpreg', 'foreach'), .errorhandling = "pass") %dopar% {
#         time.old = Sys.time()
#         fm1.j   = NULL
#         fm2.j   = NULL
#         YY.j = foreach(k = 1:K, .combine = c) %do% dat.w[[k]][,j] # scale(dat.w[[k]][,j], scale = F)
#         XX.j = matrix(0, n*K, (p-1)*K)
#         for (k in 1:K) {
#             XX.j[((k-1)*n+1):(k*n), ((k-1)*(p-1)+1):(k*(p-1))] = dat.w[[k]][,-j] # scale(dat.w[[k]][,-j], scale = F)
#         }
#         EBIC1.j = Inf
#         EBIC2.j = Inf
#         for (i.tau in 1:n.tau) {
#             if (joint) {
#                 fm.j.i = try(grpreg(XX.j, YY.j, group = group.label,
#                                     penalty = penalty, family = "gaussian",
#                                     nlambda = n.lambda, tau = taus[i.tau], lambda.min = lambda.min))
#                 if (all(class(fm.j.i) != 'try-error')) {
#                     fm.j.i = coef(fm.j.i)
#                 }
#             } else {
#                 fm.j.i = try(foreach(k = 1:K) %do% {
#                     grpreg(XX.j[((k-1)*n+1):(k*n), ((k-1)*(p-1)+1):(k*(p-1))], YY.j[((k-1)*n+1):(k*n)],
#                            group = group.label, penalty = penalty, family = "gaussian",
#                            nlambda = n.lambda, tau = taus[i.tau], lambda.min = lambda.min)
#                 })
#                 if (all(class(fm.j.i) != 'try-error')) {
#                     fm.j.i = rbind(0, foreach(k = 1:K, .combine = rbind) %do% fm.j.i[[k]]$beta[-1,])
#                 }
#             }
# 
#             if (all(class(fm.j.i) != "try-error")) {
#                 EBIC1.j.i = EBIC.wX(XX.j, YY.j, fm.j.i, probs.mat, var.common = T)
#                 EBIC2.j.i = EBIC.wX(XX.j, YY.j, fm.j.i, probs.mat, var.common = F)
#                 EBIC1.j.i[is.na(EBIC1.j.i)] = Inf
#                 EBIC2.j.i[is.na(EBIC2.j.i)] = Inf
#                 if (any(EBIC1.j.i < EBIC1.j) || any(EBIC2.j.i < EBIC2.j)) {
#                     coef.j = matrix(0, p*K, n.lambda)
#                     coef.j[-seq(j, p*K, by = p), ] = fm.j.i[-1, ]
#                     coef.j = Matrix(coef.j)
#                     if (any(EBIC1.j.i < EBIC1.j)) {
#                         fm1.j   = list(beta = coef.j, EBIC = EBIC1.j.i)
#                         EBIC1.j = min(EBIC1.j.i)
#                     }
#                     if (any(EBIC2.j.i < EBIC2.j)) {
#                         fm2.j   = list(beta = coef.j, EBIC = EBIC2.j.i)
#                         EBIC2.j = min(EBIC2.j.i)
#                     }
#                 }
#             }
#         }
#         if (verbose) { cat('.') }
#         list(fm1.j, fm2.j, Sys.time()-time.old)
#     }
# 
#     graphs.opt = foreach(crit = 1:2) %do% {
#         foreach(k = 1:K) %dopar% {
#             graph.k = foreach(j = 1:p, .combine = cbind) %do% {
#                 v.j = numeric(p)
#                 try({
#                     est.j = est.all[[j]][[crit]]
#                     if (!is.null(est.j)) {
#                         b.j.i = est.j$beta[ , which.min(est.j$EBIC)]
#                         v.j   = b.j.i[(p*(k-1)+1):(p*k)]
#                     }
#                 })
#                 v.j
#             }
#             attr(graph.k, 'dimnames') = NULL
#             Matrix(graph.k != 0)
#         }
#     }
#     graphs.all = foreach(crit = 1:2) %do% {
#         foreach(i = 1:n.lambda, .packages = 'Matrix') %dopar% {
#             foreach(k = 1:K) %do% {
#                 graph.k = foreach(j = 1:p, .combine = cbind) %do% {
#                     v.j = numeric(p)
#                     try({
#                         est.j = est.all[[j]][[crit]]
#                         if (!is.null(est.j)) {
#                             b.j.i = est.j$beta[,i]
#                             v.j   = b.j.i[(p*(k-1)+1):(p*k)]
#                         }
#                     })
#                     v.j
#                 }
#                 attr(graph.k, 'dimnames') = NULL
#                 Matrix(graph.k != 0)
#             }
#         }
#     }
#     if (verbose) { cat('\n') }
#     return(list(graphs1 = graphs.opt[[1]], graphs1.all = graphs.all[[1]],
#                 graphs2 = graphs.opt[[2]], graphs2.all = graphs.all[[2]],
#                 time = sapply(est.all, `[[`, 3), probs = probs.mat))
# }


# SoftGraphs = function(dat, label = NULL, probs.mat = NULL, joint = F, penalty = "gel",
#                       n.lambda = 20, lambda.min = 0.05, verbose = F) {
#     # Initialization
#     K       = ncol(probs.mat)
#     dat     = as.matrix(dat)
#     n       = nrow(dat)
#     p       = ncol(dat)
#     taus    = seq(0.1, 1, by = 0.1)
#     n.tau   = length(taus)
#     probs.mat = probs.mat / matrix(rowSums(probs.mat), n, K, byrow = FALSE)
#     dat.w   = ExpSample(dat, probs.mat)
# 
#     if (joint) {
#         group.label = rep(0:(p-1), K)
#     } else {
#         group.label = 0:(p-1)
#     }
# 
#     # Nodewise regression (jointly)
#     EBIC.wX = EBIC.wX
#     est.all = foreach(j = 1:p, .packages = c('grpreg', 'foreach'), .errorhandling = "pass") %dopar% {
#         time.old = Sys.time()
#         fm1.j   = NULL
#         fm2.j   = NULL
#         # try({
#             YY.j = foreach(k = 1:K, .combine = c) %do% dat.w[[k]][,j]
#             XX.j = matrix(0, n*K, p*K)
#             for (k in 1:K) {
#                 XX.j[((k-1)*n+1):(k*n), ((k-1)*p+1):(k*p)] = cbind(1, dat.w[[k]][,-j])
#             }
#             EBIC1.j = Inf
#             EBIC2.j = Inf
#             for (i.tau in 1:n.tau) {
#                 if (joint) {
#                     fm.j.i = try(grpreg(XX.j, YY.j, group = group.label,
#                                         penalty = penalty, family = "gaussian",
#                                         nlambda = n.lambda, tau = taus[i.tau], lambda.min = lambda.min))
#                     if (all(class(fm.j.i) != 'try-error')) {
#                         fm.j.i = coef(fm.j.i)
#                     }
#                 } else {
#                     fm.j.i = try(foreach(k = 1:K) %do% {
#                         grpreg(XX.j[((k-1)*n+1):(k*n), ((k-1)*p+1):(k*p)], YY.j[((k-1)*n+1):(k*n)],
#                                group = group.label, penalty = penalty, family = "gaussian",
#                                nlambda = n.lambda, tau = taus[i.tau], lambda.min = lambda.min)
#                     })
#                     if (all(class(fm.j.i) != 'try-error')) {
#                         fm.j.i = rbind(0, foreach(k = 1:K, .combine = rbind) %do% fm.j.i[[k]]$beta[-1,])
#                     }
#                 }
# 
#                 if (all(class(fm.j.i) != "try-error")) {
#                     EBIC1.j.i = EBIC.wX(XX.j, YY.j, fm.j.i, probs.mat, var.common = T)
#                     EBIC2.j.i = EBIC.wX(XX.j, YY.j, fm.j.i, probs.mat, var.common = F)
#                     EBIC1.j.i[is.na(EBIC1.j.i)] = Inf
#                     EBIC2.j.i[is.na(EBIC2.j.i)] = Inf
#                     if (any(EBIC1.j.i < EBIC1.j) || any(EBIC2.j.i < EBIC2.j)) {
#                         coef.j = matrix(0, (p+1)*K, n.lambda)
#                         coef.j[-seq(j+1, (p+1)*K, by = p+1), ] = fm.j.i[-1, ]
#                         coef.j = Matrix(coef.j[-1, ])
#                         if (any(EBIC1.j.i < EBIC1.j)) {
#                             fm1.j   = list(beta = coef.j, EBIC = EBIC1.j.i)
#                             EBIC1.j = min(EBIC1.j.i)
#                         }
#                         if (any(EBIC2.j.i < EBIC2.j)) {
#                             fm2.j   = list(beta = coef.j, EBIC = EBIC2.j.i)
#                             EBIC2.j = min(EBIC2.j.i)
#                         }
#                     }
#                 }
#             }
#         # })
#         if (verbose) { cat('.') }
#         list(fm1.j, fm2.j, Sys.time()-time.old)
#     }
# 
#     graphs.opt = foreach(crit = 1:2) %do% {
#         foreach(k = 1:K) %dopar% {
#             graph.k = foreach(j = 1:p, .combine = cbind) %do% {
#                 v.j = numeric(p)
#                 try({
#                     est.j = est.all[[j]][[crit]]
#                     if (!is.null(est.j)) {
#                         b.j.i = est.j$beta[ , which.min(est.j$EBIC)]
#                         v.j   = b.j.i[(p*(k-1)+1):(p*k)]
#                     }
#                 })
#                 v.j
#             }
#             attr(graph.k, 'dimnames') = NULL
#             Matrix(graph.k != 0)
#         }
#     }
#     graphs.all = foreach(crit = 1:2) %do% {
#         foreach(i = 1:n.lambda, .packages = 'Matrix') %dopar% {
#             foreach(k = 1:K) %do% {
#                 graph.k = foreach(j = 1:p, .combine = cbind) %do% {
#                     v.j = numeric(p)
#                     try({
#                         est.j = est.all[[j]][[crit]]
#                         if (!is.null(est.j)) {
#                             b.j.i = est.j$beta[,i]
#                             v.j   = b.j.i[(p*(k-1)+1):(p*k)]
#                         }
#                     })
#                     v.j
#                 }
#                 attr(graph.k, 'dimnames') = NULL
#                 Matrix(graph.k != 0)
#             }
#         }
#     }
#     if (verbose) { cat('\n') }
#     return(list(graphs1 = graphs.opt[[1]], graphs1.all = graphs.all[[1]],
#                 graphs2 = graphs.opt[[2]], graphs2.all = graphs.all[[2]],
#                 time = sapply(est.all, `[[`, 3), probs = probs.mat))
# }

SoftGraphs = function(dat, label = NULL, probs.mat = NULL, joint = F, penalty = "gel",
                      n.lambda = 20, lambda.min = 0.05, verbose = F) {
    if (is.null(probs.mat)) {
        stopifnot(!is.null(label))
        probs.mat = SoftLabel.QDA(dat, label)
    }
    K.orig = ncol(probs.mat)
    if (all(probs.mat == probs.mat[,1])) {
        probs.mat = matrix(1, nrow(probs.mat), 1)
    }
    # Initialization
    K       = ncol(probs.mat)
    dat     = as.matrix(dat)
    n       = nrow(dat)
    p       = ncol(dat)
    taus    = seq(0.1, 1, by = 0.1)
    n.tau   = length(taus)
    probs.mat = probs.mat / matrix(rowSums(probs.mat), n, K, byrow = FALSE)
    probs.mat[probs.mat < 1e-6] = 0
    dat.w   = ExpSample(dat, probs.mat, thres = 1e-6)
    n.ext   = sum(sapply(dat.w, nrow))
    obs.grp = cumsum(sapply(dat.w, nrow))
    obs.grp = cbind(c(1, obs.grp[-K]+1), obs.grp)
    obs.grp = foreach(k = 1:K) %do% obs.grp[k,1]:obs.grp[k,2]

    if (joint) {
        group.label = rep(1:(p-1), K)
    } else {
        group.label = 1:(p-1)
    }

    # Nodewise regression (jointly)
    EBIC.wX = EBIC.wX
    est.all = foreach(j = 1:p, .packages = c('grpreg', 'foreach'), .errorhandling = "pass") %dopar% {
      library(Matrix)
        time.old = Sys.time()
        fm1.j   = NULL
        fm2.j   = NULL
        try({
        YY.j = foreach(k = 1:K, .combine = c) %do% scale(dat.w[[k]][,j], scale = F)
        XX.j = matrix(0, n.ext, K*(p-1))
        for (k in 1:K) {
            XX.j[obs.grp[[k]], ((k-1)*(p-1)+1):(k*(p-1))] = scale(dat.w[[k]][,-j], scale = F)
        }
        EBIC1.j = Inf
        EBIC2.j = Inf
        for (i.tau in 1:n.tau) {
            if (joint) {
                fm.j.i = try(grpreg(XX.j, YY.j, group = group.label,
                                    penalty = penalty, family = "gaussian",
                                    nlambda = n.lambda, tau = taus[i.tau], lambda.min = lambda.min))
                if (all(class(fm.j.i) != 'try-error')) {
                    fm.j.i = coef(fm.j.i)
                }
            } else {
                fm.j.i = try(foreach(k = 1:K) %do% {
                    grpreg(XX.j[obs.grp[[k]], ((k-1)*(p-1)+1):(k*(p-1))], YY.j[obs.grp[[k]]],
                           group = group.label, penalty = penalty, family = "gaussian",
                           nlambda = n.lambda, tau = taus[i.tau], lambda.min = lambda.min)
                })
                if (all(class(fm.j.i) != 'try-error')) {
                    fm.j.i = rbind(0, foreach(k = 1:K, .combine = rbind) %do% fm.j.i[[k]]$beta[-1,])
                }
            }

            if (all(class(fm.j.i) != "try-error")) {
                EBIC1.j.i = EBIC.wX(XX.j, YY.j, fm.j.i, probs.mat, var.common = T)
                EBIC2.j.i = EBIC.wX(XX.j, YY.j, fm.j.i, probs.mat, var.common = F, group = obs.grp)
                EBIC1.j.i[is.na(EBIC1.j.i)] = Inf
                EBIC2.j.i[is.na(EBIC2.j.i)] = Inf
                if (any(EBIC1.j.i < EBIC1.j) || any(EBIC2.j.i < EBIC2.j)) {
                    coef.j = matrix(0, p*K, n.lambda)
                    coef.j[-seq(j, p*K, by = p), ] = fm.j.i[-1, ]
                    coef.j = Matrix(coef.j)
                    if (any(EBIC1.j.i < EBIC1.j)) {
                        fm1.j   = list(beta = coef.j, EBIC = EBIC1.j.i)
                        EBIC1.j = min(EBIC1.j.i)
                    }
                    if (any(EBIC2.j.i < EBIC2.j)) {
                        fm2.j   = list(beta = coef.j, EBIC = EBIC2.j.i)
                        EBIC2.j = min(EBIC2.j.i)
                    }
                }
            }
        }
        })
        if (verbose) { cat('.') }
        list(fm1.j, fm2.j, Sys.time()-time.old)
    }

    graphs.opt = foreach(crit = 1:2) %do% {
        foreach(k = 1:K) %dopar% {
          library(Matrix)
            graph.k = foreach(j = 1:p, .combine = cbind) %do% {
                v.j = numeric(p)
                try({
                    est.j = est.all[[j]][[crit]]
                    if (!is.null(est.j)) {
                        b.j.i = est.j$beta[ , which.min(est.j$EBIC)]
                        v.j   = b.j.i[(p*(k-1)+1):(p*k)]
                    }
                })
                v.j
            }
            attr(graph.k, 'dimnames') = NULL
            Matrix(graph.k != 0)
        }
    }
    graphs.all = foreach(crit = 1:2) %do% {
        foreach(i = 1:n.lambda, .packages = 'Matrix') %dopar% {
            foreach(k = 1:K) %do% {
              library(Matrix)
                graph.k = foreach(j = 1:p, .combine = cbind) %do% {
                    v.j = numeric(p)
                    try({
                        est.j = est.all[[j]][[crit]]
                        if (!is.null(est.j)) {
                            b.j.i = est.j$beta[,i]
                            v.j   = b.j.i[(p*(k-1)+1):(p*k)]
                        }
                    })
                    v.j
                }
                attr(graph.k, 'dimnames') = NULL
                Matrix(graph.k != 0)
            }
        }
    }
    if (K.orig > K) {
        graphs.opt[[1]] = rep(graphs.opt[[1]], K.orig)
        graphs.opt[[2]] = rep(graphs.opt[[2]], K.orig)
        graphs.all[[1]] = lapply(graphs.all[[1]], function(x) rep(x, K.orig))
        graphs.all[[2]] = lapply(graphs.all[[2]], function(x) rep(x, K.orig))
        probs.mat = matrix(1/K.orig, n, K.orig)
    }

    if (verbose) { cat('\n') }
    return(list(graphs1 = graphs.opt[[1]], graphs1.all = graphs.all[[1]],
                graphs2 = graphs.opt[[2]], graphs2.all = graphs.all[[2]],
                time = sapply(est.all, `[[`, 3), probs = probs.mat))
}


SymmetrizeSoftGraphs = function(file.res, method = c('Union', 'Intersect')) {
    method = match.arg(method)
    newfile.res = paste0('Symm', method, '_', file.res)
    
    load(file.res)
    try({res.st = SymmetrizeRes(res.st, method)})
    try({res.hs = SymmetrizeRes(res.hs, method)})
    try({res.ss = SymmetrizeRes(res.ss, method)})
    try({res.hj = SymmetrizeRes(res.hj, method)})
    try({res.sj = SymmetrizeRes(res.sj, method)})
    
    save(train, K, genes.now, res.st, res.hs, res.ss, res.hj, res.sj, genefile, 
         file = newfile.res)
    cat('Graph symmetrized by', method, '. Saved to ', newfile.res, '\n')
    gc(verbose = FALSE)
    return(0)
}

SymmetrizeRes = function(res, method = c('Union', 'Intersect')) {
    method = match.arg(method)
    tmp = list()
    if (method == 'Union') {
        tmp$graphs1 = lapply(res$graphs1, 
            function(x) {x = as.matrix(x); attr(x, 'dimnames') = NULL; Matrix(x | t(x))})
        tmp$graphs2 = lapply(res$graphs2, 
            function(x) {x = as.matrix(x); attr(x, 'dimnames') = NULL; Matrix(x | t(x))})
        tmp$graphs1.all = lapply(res$graphs1.all, 
            function(xs) lapply(xs, 
                function(x) {x = as.matrix(x); attr(x, 'dimnames') = NULL; Matrix(x | t(x))}))
        tmp$graphs2.all = lapply(res$graphs2.all, 
            function(xs) lapply(xs, 
                function(x) {x = as.matrix(x); attr(x, 'dimnames') = NULL; Matrix(x | t(x))}))
    } else {
        tmp$graphs1 = lapply(res$graphs1, 
                             function(x) {x = as.matrix(x); attr(x, 'dimnames') = NULL; Matrix(x & t(x))})
        tmp$graphs2 = lapply(res$graphs2, 
                             function(x) {x = as.matrix(x); attr(x, 'dimnames') = NULL; Matrix(x & t(x))})
        tmp$graphs1.all = lapply(res$graphs1.all, 
            function(xs) lapply(xs, 
                function(x) {x = as.matrix(x); attr(x, 'dimnames') = NULL; Matrix(x & t(x))}))
        tmp$graphs2.all = lapply(res$graphs2.all, 
            function(xs) lapply(xs, 
                function(x) {x = as.matrix(x); attr(x, 'dimnames') = NULL; Matrix(x & t(x))}))
    }
    tmp$probs = res$probs
    return(tmp)
}

# Weighted JGL
# Y is a matrix with each row for an observation
WJGL = function(Y, probs, penalty = "fused", lambda1, lambda2, 
                rho = 1, weights = "equal", penalize.diagonal = FALSE, 
                maxiter = 500, tol = 1e-05, warm = NULL, 
                return.whole.theta = FALSE, screening = "fast", 
                truncate = 1e-05) {
    
    K     = ncol(probs)
    
    dat.w = WSample(Y, probs)
    
    if (length(weights) == 1) {
        stopifnot(weights == "equal")
        weights = dat.w$weights
    } else {
        stopifnot(length(weights) == K)
        weights = dat.w$weights * weights
        weights = weights / sum(weights)
    }
    
    output = JGL(dat.w$Y, penalty, lambda1, lambda2, rho, weights, 
                 penalize.diagonal, maxiter, tol, warm, 
                 return.whole.theta, screening, truncate)
    
    return(output)
}
