
# PC algorithm, the stable version.
pc.stable = function(x, cluster, whitelist, blacklist, test = NULL,
    alpha = 0.05, ..., max.sx = NULL, debug = FALSE, undirected = FALSE) {

  bnlearn(data = x, cluster = cluster, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha, extra.args = list(...),
    max.sx = max.sx, algorithm = "pc.stable", debug = debug,
    undirected = undirected)

}#PC.STABLE

# Grow-Shrink frontend.
gs = function(x, cluster, whitelist, blacklist, test = NULL, alpha = 0.05, ...,
    max.sx = NULL, debug = FALSE, undirected = FALSE) {

  bnlearn(data = x, cluster = cluster, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha, extra.args = list(...),
    max.sx = max.sx, algorithm = "gs", debug = debug, undirected = undirected)

}#GS

# Incremental Association frontend.
iamb = function(x, cluster, whitelist, blacklist, test = NULL, alpha = 0.05,
    ..., max.sx = NULL, debug = FALSE, undirected = FALSE) {

  bnlearn(data = x, cluster = cluster, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha, extra.args = list(...),
    max.sx = max.sx, algorithm = "iamb", debug = debug, undirected = undirected)

}#IAMB

# Fast-IAMB frontend.
fast.iamb = function(x, cluster, whitelist, blacklist, test = NULL,
    alpha = 0.05, ..., max.sx = NULL, debug = FALSE, undirected = FALSE) {

  bnlearn(data = x, cluster = cluster, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha, extra.args = list(...),
    max.sx = max.sx, algorithm = "fast.iamb", debug = debug,
    undirected = undirected)

}#FAST.IAMB

# Inter-IAMB frontend.
inter.iamb = function(x, cluster, whitelist, blacklist, test = NULL,
    alpha = 0.05, ...,  max.sx = NULL, debug = FALSE, undirected = FALSE) {

  bnlearn(data = x, cluster = cluster, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha, extra.args = list(...),
    max.sx = max.sx, algorithm = "inter.iamb", debug = debug,
    undirected = undirected)

}#INTER.IAMB

# IAMB-FDR frontend.
iamb.fdr = function(x, cluster, whitelist, blacklist, test = NULL, alpha = 0.05,
    ...,  max.sx = NULL, debug = FALSE, undirected = FALSE) {

  bnlearn(data = x, cluster = cluster, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha, extra.args = list(...),
    max.sx = max.sx, algorithm = "iamb.fdr", debug = debug,
    undirected = undirected)

}#IAMB.FDR

# MMPC frontend.
mmpc = function(x, cluster, whitelist, blacklist, test = NULL, alpha = 0.05,
    ..., max.sx = NULL, debug = FALSE, undirected = TRUE) {

  bnlearn(data = x, cluster = cluster, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha, extra.args = list(...),
    max.sx = max.sx, algorithm = "mmpc", debug = debug,
    undirected = undirected)

}#MMPC

# Semi-Interleaved HITON-PC.
si.hiton.pc = function(x, cluster, whitelist, blacklist, test = NULL,
    alpha = 0.05, ..., max.sx = NULL, debug = FALSE, undirected = TRUE) {

  bnlearn(data = x, cluster = cluster, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha, extra.args = list(...),
    max.sx = max.sx, algorithm = "si.hiton.pc", debug = debug,
    undirected = undirected)

}#SI.HITON.PC

# Hybrid PC.
hpc = function(x, cluster, whitelist, blacklist, test = NULL, alpha = 0.05,
    ..., max.sx = NULL, debug = FALSE, undirected = TRUE) {

  bnlearn(data = x, cluster = cluster, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha, extra.args = list(...),
    max.sx = max.sx, algorithm = "hpc", debug = debug, undirected = undirected)

}#HPC

# ARACNE frontend.
aracne = function(x, whitelist, blacklist, mi = NULL, debug = FALSE) {

  mi.matrix(data = x, whitelist = whitelist, blacklist = blacklist,
    algorithm = "aracne", mi = mi, debug = debug)

}#ARACNE

# Chow-Liu frontend.
chow.liu  = function(x, whitelist, blacklist, mi = NULL, debug = FALSE) {

  mi.matrix(data = x, whitelist = whitelist, blacklist = blacklist,
    algorithm = "chow.liu", mi = mi, debug = debug)

}#CHOW.LIU

# Hill Climbing greedy search frontend.
hc = function(x, start = NULL, whitelist, blacklist, score = NULL, ...,
    debug = FALSE, restart = 0, perturb = 1, max.iter = Inf, maxp = Inf,
    optimized = TRUE) {

  greedy.search(data = x, start = start, whitelist = whitelist,
    blacklist = blacklist, score = score, algorithm = "hc", debug = debug,
    ..., restart = restart, perturb = perturb,
    max.iter = max.iter, maxp = maxp, optimized = optimized)

}#HC

# TABU list greedy search frontend.
tabu = function(x, start = NULL, whitelist, blacklist, score = NULL, ...,
    debug = FALSE, tabu = 10, max.tabu = tabu, max.iter = Inf, maxp = Inf,
    optimized = TRUE) {

  greedy.search(data = x, start = start, whitelist = whitelist,
    blacklist = blacklist, score = score, algorithm = "tabu", debug = debug,
    ..., max.iter = max.iter, tabu = tabu, max.tabu = max.tabu,
    maxp = maxp, optimized = optimized)

}#TABU

# Generic Restricted Maximization frontend.
rsmax2 = function(x, whitelist, blacklist, restrict = "si.hiton.pc",
    maximize = "hc", restrict.args = list(), maximize.args = list(),
    debug = FALSE) {

  hybrid.search(data = x, whitelist = whitelist, blacklist = blacklist,
    restrict = restrict, maximize = maximize, restrict.args = restrict.args,
    maximize.args = maximize.args, debug = debug)

}#RSMAX2

# MMHC frontend.
mmhc = function(x, whitelist, blacklist, restrict.args = list(),
    maximize.args = list(), debug = FALSE) {

  hybrid.search(data = x, whitelist = whitelist, blacklist = blacklist,
    restrict = "mmpc", maximize = "hc", restrict.args = restrict.args,
    maximize.args = maximize.args, debug = debug)

}#MMHC

# H2PC frontend.
h2pc = function(x, whitelist, blacklist, restrict.args = list(),
    maximize.args = list(), debug = FALSE) {

  hybrid.search(data = x, whitelist = whitelist, blacklist = blacklist,
    restrict = "hpc", maximize = "hc", restrict.args = restrict.args,
    maximize.args = maximize.args, debug = debug)

}#H2PC

# Frontend for the Markov blanket learning algorithms.
learn.mb = function(x, node, method, whitelist, blacklist, start = NULL,
    test = NULL, alpha = 0.05, ..., max.sx = NULL, debug = FALSE) {

  mb.backend(data = x, target = node, algorithm = method, whitelist = whitelist,
    blacklist = blacklist, start = start, test = test, alpha = alpha,
    extra.args = list(...), max.sx = max.sx, debug = debug)

}#LEARN.MB

# Frontend for causal discovery learning algorithms.
learn.nbr = function(x, node, method, whitelist, blacklist, test = NULL,
    alpha = 0.05, ..., max.sx = NULL, debug = FALSE) {

  nbr.backend(data = x, target = node, algorithm = method,
    whitelist = whitelist, blacklist = blacklist, test = test, alpha = alpha,
    extra.args = list(...), max.sx = max.sx, debug = debug)

}#LEARN.NBR

# naive Bayes frontend.
naive.bayes = function(x, training, explanatory) {

  bayesian.classifier(data = x, training = training, explanatory = explanatory,
    algorithm = "naive.bayes", whitelist = NULL, blacklist = NULL,
    expand = list(), debug = FALSE)

}#NAIVE.BAYES

# tree-augmented naive Bayes frontend.
tree.bayes = function(x, training, explanatory, whitelist, blacklist,
    mi = NULL, root = NULL, debug = FALSE) {

  bayesian.classifier(data = x, training = training, explanatory = explanatory,
    algorithm = "tree.bayes", whitelist = whitelist, blacklist = blacklist,
    expand = list(estimator = mi, root = root), debug = debug)

}#TREE.BAYES

# Direct LiNGAM frontend.
direct.lingam = function(x, cluster, whitelist, blacklist, mi,
    maximize = "alasso", maximize.args = list(), debug = FALSE) {

  lingam.learners(data = x, cluster = cluster, algorithm = "direct.lingam",
    whitelist = whitelist, blacklist = blacklist, mi = mi, maximize = maximize,
    maximize.args = maximize.args, debug = debug)

}#DIRECT.LINGAM
