
# PC algorithm, the stable version.
pc.stable = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
    test = NULL, alpha = 0.05, B = NULL, max.sx = NULL, debug = FALSE,
    undirected = FALSE) {

  bnlearn(x = x, cluster = cluster, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha, B = B,
    max.sx = max.sx, method = "pc.stable", debug = debug,
    undirected = undirected)

}#PC.CLASSIC

# Grow-Shrink frontend.
gs = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
    test = NULL, alpha = 0.05, B = NULL, max.sx = NULL, debug = FALSE,
    undirected = FALSE) {

  bnlearn(x = x, cluster = cluster, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha, B = B,
    max.sx = max.sx, method = "gs", debug = debug, undirected = undirected)

}#GS

# Incremental Association frontend.
iamb = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
    test = NULL, alpha = 0.05, B = NULL, max.sx = NULL, debug = FALSE,
    undirected = FALSE) {

  bnlearn(x = x, cluster = cluster, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha, B = B,
    max.sx = max.sx, method = "iamb", debug = debug, undirected = undirected)

}#IAMB

# Fast-IAMB frontend.
fast.iamb = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
    test = NULL, alpha = 0.05, B = NULL, max.sx = NULL, debug = FALSE,
    undirected = FALSE) {

  bnlearn(x = x, cluster = cluster, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha, B = B,
    max.sx = max.sx, method = "fast.iamb", debug = debug,
    undirected = undirected)

}#FAST.IAMB

# Inter-IAMB frontend.
inter.iamb = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
    test = NULL, alpha = 0.05, B = NULL,  max.sx = NULL, debug = FALSE,
    undirected = FALSE) {

  bnlearn(x = x, cluster = cluster, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha, B = B,
    max.sx = max.sx, method = "inter.iamb", debug = debug,
    undirected = undirected)

}#INTER.IAMB

# IAMB-FDR frontend.
iamb.fdr = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
    test = NULL, alpha = 0.05, B = NULL,  max.sx = NULL, debug = FALSE,
    undirected = FALSE) {

  bnlearn(x = x, cluster = cluster, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha, B = B,
    max.sx = max.sx, method = "iamb.fdr", debug = debug,
    undirected = undirected)

}#IAMB.FDR

# MMPC frontend.
mmpc = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
    test = NULL, alpha = 0.05, B = NULL, max.sx = NULL, debug = FALSE,
    undirected = TRUE) {

  bnlearn(x = x, cluster = cluster, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha, B = B,
    max.sx = max.sx, method = "mmpc", debug = debug,
    undirected = undirected)

}#MMPC

# Semi-Interleaved HITON-PC.
si.hiton.pc = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
    test = NULL, alpha = 0.05, B = NULL, max.sx = NULL, debug = FALSE,
    undirected = TRUE) {

  bnlearn(x = x, cluster = cluster, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha, B = B,
    max.sx = max.sx, method = "si.hiton.pc", debug = debug,
    undirected = undirected)

}#SI.HITON.PC

# Hybrid PC.
hpc = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
    test = NULL, alpha = 0.05, B = NULL, max.sx = NULL, debug = FALSE,
    undirected = TRUE) {

  bnlearn(x = x, cluster = cluster, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha, B = B,
    max.sx = max.sx, method = "hpc", debug = debug, undirected = undirected)

}#HPC

# ARACNE frontend.
aracne = function(x, whitelist = NULL, blacklist = NULL, mi = NULL,
    debug = FALSE) {

  mi.matrix(x = x, whitelist = whitelist, blacklist = blacklist,
    method = "aracne", mi = mi, debug = debug)

}#ARACNE

# Chow-Liu frontend.
chow.liu  = function(x, whitelist = NULL, blacklist = NULL, mi = NULL,
    debug = FALSE) {

  mi.matrix(x = x, whitelist = whitelist, blacklist = blacklist,
    method = "chow.liu", mi = mi, debug = debug)

}#CHOW.LIU

# Hill Climbing greedy search frontend.
hc = function(x, start = NULL, whitelist = NULL, blacklist = NULL,
    score = NULL, ..., debug = FALSE, restart = 0, perturb = 1,
    max.iter = Inf, maxp = Inf, optimized = TRUE) {

  greedy.search(x = x, start = start, whitelist = whitelist,
    blacklist = blacklist, score = score, heuristic = "hc", debug = debug,
    ..., restart = restart, perturb = perturb,
    max.iter = max.iter, maxp = maxp, optimized = optimized)

}#HC

# TABU list greedy search frontend.
tabu = function(x, start = NULL, whitelist = NULL, blacklist = NULL,
    score = NULL, ..., debug = FALSE, tabu = 10, max.tabu = tabu,
    max.iter = Inf, maxp = Inf, optimized = TRUE) {

  greedy.search(x = x, start = start, whitelist = whitelist,
    blacklist = blacklist, score = score, heuristic = "tabu", debug = debug,
    ..., max.iter = max.iter, tabu = tabu, max.tabu = max.tabu,
    maxp = maxp, optimized = optimized)

}#TABU

# Generic Restricted Maximization frontend.
rsmax2 = function(x, whitelist = NULL, blacklist = NULL, restrict = "si.hiton.pc",
    maximize = "hc", restrict.args = list(), maximize.args = list(),
    debug = FALSE) {

  hybrid.search(x, whitelist = whitelist, blacklist = blacklist,
    restrict = restrict, maximize = maximize, restrict.args = restrict.args,
    maximize.args = maximize.args, debug = debug)

}#RSMAX2

# MMHC frontend.
mmhc = function(x, whitelist = NULL, blacklist = NULL,
    restrict.args = list(), maximize.args = list(), debug = FALSE) {

  hybrid.search(x, whitelist = whitelist, blacklist = blacklist,
    restrict = "mmpc", maximize = "hc", restrict.args = restrict.args,
    maximize.args = maximize.args, debug = debug)

}#MMHC

# H2PC frontend.
h2pc = function(x, whitelist = NULL, blacklist = NULL,
    restrict.args = list(), maximize.args = list(), debug = FALSE) {

  hybrid.search(x, whitelist = whitelist, blacklist = blacklist,
    restrict = "hpc", maximize = "hc", restrict.args = restrict.args,
    maximize.args = maximize.args, debug = debug)

}#H2PC

# Frontend for the Markov blanket learning algorithms.
learn.mb = function(x, node, method, whitelist = NULL, blacklist = NULL,
    start = NULL, test = NULL, alpha = 0.05, B = NULL, max.sx = NULL,
    debug = FALSE) {

  mb.backend(x, target = node, method = method, whitelist = whitelist,
    blacklist = blacklist, start = start, test = test, alpha = alpha,
    B = B, max.sx = max.sx, debug = debug)

}#LEARN.MB

# Frontend for causal discovery learning algorithms.
learn.nbr = function(x, node, method, whitelist = NULL, blacklist = NULL,
    test = NULL, alpha = 0.05, B = NULL, max.sx = NULL, debug = FALSE) {

  nbr.backend(x, target = node, method = method, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha, B = B,
    max.sx = max.sx, debug = debug)

}#LEARN.NBR

# naive Bayes frontend.
naive.bayes = function(x, training, explanatory) {

  bayesian.classifier(x, training = training, explanatory = explanatory,
    method = "naive.bayes", whitelist = NULL, blacklist = NULL, expand = list(),
    debug = FALSE)

}#NAIVE.BAYES

# tree-augmented naive Bayes frontend.
tree.bayes = function(x, training, explanatory, whitelist = NULL, blacklist = NULL,
    mi = NULL, root = NULL, debug = FALSE) {

  bayesian.classifier(x, training = training, explanatory = explanatory,
    method = "tree.bayes", whitelist = whitelist, blacklist = blacklist,
    expand = list(estimator = mi, root = root), debug = debug)

}#TREE.BAYES
