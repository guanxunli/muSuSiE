# Load Packages and Programs
source("JDAG_head.r")
library(HDclassif)

# # Parallel Computing
# cl = makeCluster(4)
# registerDoParallel(cl)
# getDoParWorkers()

# Settings
graph.type     = 'BA'           # BA or ER
K              = 4              # no. of classes (K)
p              = 500            # dimension (p)
n              = 100            # average sample size
ndelta         = 4              
sps            = 1              # graph sparsity
p.common       = 0.7            # graph overlapping (\pi_0)
diff.scale     = 0.05           # mean difference level (\delta^2)
n.graph        = 50             # no. of estimated graphs
lambda.min     = 0.05           # lower bound of regularization parameter

set.seed(41)

# Precision Matrices and Graphs
OmegaDAG   = SimOmegas(K, p, type = graph.type, sparsity = sps/p, p.common = p.common, 
                       common.strength = 1, unique.strength = 1, 
                       u.seed = sample(1e8, 1), disp = F)
Gammas     = OmegaDAG$DAGs
Omegas     = OmegaDAG$Omegas
DAGs       = lapply(Gammas, function(x) {diag(x) = 0; 0 + (t(x) != 0)})
CPDAGs.g   = lapply(DAGs, function(x) graph_from_adjacency_matrix(x, mode = 'directed', diag=F))
CPDAGs.g   = lapply(CPDAGs.g, function(x) dag2cpdag(as_graphnel(x)))
CPDAGs     = lapply(CPDAGs.g, graph2Matrix, symmetric = F)
skeletons  = lapply(DAGs, function(x) {y = AsSymmetric(x); diag(y) = 0; 0 + (y != 0)})

# Means
p.uni       = ndelta
mu.all      = matrix(0, K, p)
for (k in 1:K) {
    mu.all[k, ((k-1)*p.uni + 1):(k*p.uni)] = sqrt(diff.scale)
}
ix.uni      = 1:(K*p.uni)
nu.all      = foreach(k = 1:K, .combine = rbind) %do% {
    c(Gammas[[k]] %*% mu.all[k, ])
}
n.train    = sample(K, n*K, replace=T)
n.train <- as.numeric(table(n.train))
train      = SimGaussian.B(n.train, nu.all, Gammas)
N.train <- sum(n.train)
probs.act  = ProbsQDA(train$X, mu.all, Omegas, w = n.train / N.train)
colnames(probs.act) = as.character(1:K)

# Label Generating
svd.X = svd(train$X)
pc.X = svd.X$u %*% diag(svd.X$d)[,1:min(20, p)]
# Avoid Unreasonable Clustering
for (i in 1:100) {
    res.cl = RegisterLabel(train$g, hddc(pc.X, K = K)$class)
             # kmeans(dist(train$X[,ix.sig]), K)$cluster)
    cat('Error: ', round(mean(res.cl != train$g), 4))
    if (mean(res.cl != train$g) < 0.2) { break }
}

hlab.cl = matrix(0, N.train, K)
hlab.cl[cbind(1:N.train, res.cl)] = 1
hlab.oracle = matrix(0, N.train, K)
hlab.oracle[cbind(1:N.train, train$g)] = 1

slab.oracle = probs.act
slab.NB     = SoftLabel.QDA(pc.X, res.cl)
slab.triv   = matrix(1 / K, nrow = N.train, ncol = K)


# Stage I.  Undirected Graph Estimation
# PenPC w/o Grouping
res.st <- try(SoftGraphs(train$X, probs.mat = slab.triv, joint = F, 
                         n.lambda = n.graph, lambda.min = lambda.min))
# Hard PenPC
res.hs <- try(SoftGraphs(train$X, probs.mat = hlab.cl, joint = F, 
                         n.lambda = n.graph, lambda.min = lambda.min))
# Soft PenPC
res.ss <- try(SoftGraphs(train$X, probs.mat = slab.NB, joint = F, 
                         n.lambda = n.graph, lambda.min = lambda.min))
# Hard MPenPC
res.hj <- try(SoftGraphs(train$X, probs.mat = hlab.cl, joint = T, 
                         n.lambda = n.graph, lambda.min = lambda.min))
# Soft MPenPC
res.sj <- try(SoftGraphs(train$X, probs.mat = slab.NB, joint = T, 
                         n.lambda = n.graph, lambda.min = lambda.min))
PEN.list = list(res.st, res.hs, res.ss, res.hj, res.sj)


# Stage II. Skeleton Estimation (Fixed alpha)
PEN.list2 = lapply(PEN.list, SymmetrizeRes, 'Union')
ini.list  = lapply(PEN.list2, function(x) list(graphs = x$graphs1, probs = x$probs))
PC.list   = foreach(i = 1:length(ini.list), .packages = c("foreach", "pcalg")) %dopar% {
    library(ParallelPC, lib = c(.libPaths(), "/nas02/home/l/i/liuoo/Rlibs"))
    source("JDAG_Functions_Test.r")
    time.old = Sys.time()
    tmp = SoftPC(ini.list[[i]], X = train$X, alpha = 0.02)
    tmp$time = Sys.time() - time.old
    tmp
}


# Stage II. Skeleton Estimation (Varying alpha)
alphas = seq(0.005, 0.03, by = 0.005)
PC.lists = foreach(alpha = alphas) %do% {
    foreach(i = 1:length(ini.list), .packages = c("foreach", "pcalg")) %dopar% {
        library(ParallelPC, lib = c(.libPaths(), "/nas02/home/l/i/liuoo/Rlibs"))
        source("JDAG_Functions_Test.r")
        time.old = Sys.time()
        tmp = SoftPC(ini.list[[i]], X = train$X, alpha = alpha)
        tmp$time = Sys.time() - time.old
        tmp
    }
}


# Evaluation of Selected Estimation
perf.PEN = foreach(i = 1:length(PEN.list), .combine = rbind, .packages = 'foreach') %dopar% {
    MeanGraphsAccuracy(PEN.list[[i]]$graphs1, skeletons)
}

perf.PC = foreach(i = 1:length(PC.list), .combine = rbind) %do% {
    PCs.i = foreach(k = 1:K) %do% graph2Matrix(PC.list[[i]][[k]]@graph)
    MeanGraphsAccuracy(PCs.i, skeletons)
}

perf.shd = foreach(i = 1:length(PC.list), .combine = c) %do% {
    cpdags.i = lapply(PC.list[[i]][1:K], udag2pdag)
    mean(mapply(shd, cpdags.i, CPDAGs.g))
}
perf.shd = data.frame(shd = perf.shd, Method = mths.disp, pi0 = p.common, piE = sps, 
                      b = str.Omega, delta = diff.scale, n = n, nd = ndelta)

# Stage I Evaluation
perf.vPEN = foreach(v = 1:length(PEN.list), .combine = rbind, .packages = 'foreach') %dopar% {
    foreach(i = 1:n.graph, .combine = rbind) %do% {
        MeanGraphsAccuracy(PEN.list[[v]]$graphs1.all[[i]], skeletons)
    }
}

# Stage II Evaluation
perf.vPC = foreach(i = 1:length(ini.list), .packages = 'foreach') %do% {
    perf.i = foreach(i.alpha = 1:length(PC.lists), .combine = rbind) %do% {
        skeletons.ii = lapply(PC.lists[[i.alpha]][[i]][1:K], function(x) graph2Matrix(x@graph))
        cbind(MeanGraphsAccuracy(skeletons.ii, skeletons), time = PC.lists[[i.alpha]][[i]]$time)
    }
}







