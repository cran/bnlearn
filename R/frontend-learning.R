
# Grow-Shrink frontend.
gs = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
    test = NULL, alpha = 0.05, B = NULL, debug = FALSE, optimized = TRUE,
    strict = FALSE, undirected = FALSE) {

  bnlearn(x = x, cluster = cluster, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha, B = B, debug = debug,
    optimized = optimized, strict = strict, undirected = undirected)

}#GS

# Incremental Association frontend.
iamb = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
    test = NULL, alpha = 0.05, B = NULL, debug = FALSE, optimized = TRUE,
    strict = FALSE, undirected = FALSE) {

  bnlearn(x = x, cluster = cluster, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha, B = B, method = "iamb",
    debug = debug, optimized = optimized, strict = strict,
    undirected = undirected)

}#IAMB

# Fast-IAMB frontend.
fast.iamb = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
    test = NULL, alpha = 0.05, B = NULL, debug = FALSE, optimized = TRUE,
    strict = FALSE, undirected = FALSE) {

  bnlearn(x = x, cluster = cluster, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha, B = B,
    method = "fast.iamb", debug = debug, optimized = optimized,
    strict = strict, undirected = undirected)

}#FAST.IAMB

# Inter-IAMB frontend.
inter.iamb = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
    test = NULL, alpha = 0.05, B = NULL, debug = FALSE, optimized = TRUE,
    strict = FALSE, undirected = FALSE) {

  bnlearn(x = x, cluster = cluster, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha, B = B,
    method = "inter.iamb", debug = debug, optimized = optimized,
    strict = strict, undirected =  undirected)

}#INTER.IAMB

# MMPC frontend.
mmpc = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
    test = NULL, alpha = 0.05, B = NULL, debug = FALSE, optimized = TRUE,
    strict = FALSE) {

  bnlearn(x = x, cluster = cluster, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha, B = B,
    method = "mmpc", debug = debug, optimized = optimized,
    strict = strict, undirected = TRUE)

}#MMPC

# Hill Climbing greedy search frontend.
hc = function(x, start = NULL, whitelist = NULL, blacklist = NULL,
    score = NULL, ..., debug = FALSE, restart = 0, perturb = 1,
    max.iter = Inf, optimized = TRUE) {

  greedy.search(x = x, start = start, whitelist = whitelist,
    blacklist = blacklist, score = score, heuristic = "hc", ...,
    debug = debug, misc.args = list(restart = restart, perturb = perturb,
    max.iter = max.iter), optimized = optimized)

}#HC

# TABU list greedy search frontend.
tabu = function(x, start = NULL, whitelist = NULL, blacklist = NULL,
    score = NULL, ..., debug = FALSE, tabu = 10, max.tabu = tabu,
    max.iter = Inf, optimized = TRUE) {

  greedy.search(x = x, start = start, whitelist = whitelist,
    blacklist = blacklist, score = score, heuristic = "tabu", ...,
    debug = debug, misc.args = list(max.iter = max.iter, tabu = tabu,
    max.tabu = max.tabu), optimized = optimized)

}#TABU

# Generic Restricted Maximization frontend.
rsmax2 = function(x, whitelist = NULL, blacklist = NULL, restrict = "gs",
    maximize = "hc", test = NULL, score = NULL, alpha = 0.05, B = NULL,
    ..., maximize.args = list(), optimized = TRUE, strict = FALSE,
    debug = FALSE) {

  hybrid.search(x, whitelist = whitelist, blacklist = blacklist,
    restrict = restrict, maximize = maximize, test = test, score = score,
    alpha = alpha, B = B, ..., maximize.args = maximize.args,
    optimized = optimized, strict = strict, debug = debug)

}#RSHC

# MMHC frontend.
mmhc = function(x, whitelist = NULL, blacklist = NULL, test = NULL,
    score = NULL, alpha = 0.05, B = NULL, ..., restart = 0, perturb = 1,
    max.iter = Inf, optimized = TRUE, strict = FALSE, debug = FALSE) {

  max.args = list(restart = restart, perturb = perturb, max.iter = max.iter)

  hybrid.search(x, whitelist = whitelist, blacklist = blacklist,
    restrict = "mmpc", maximize = "hc", test = test, score = score,
    alpha = alpha, B = B, ..., maximize.args = max.args,
    optimized = optimized, strict = strict,
    debug = debug)

}#MMHC

# Frontend for the Markov blanket learning algotrithms.
learn.mb = function(x, node, method, whitelist = NULL, blacklist = NULL,
    test = NULL, alpha = 0.05, B = NULL, debug = FALSE, optimized = TRUE) {

  mb.backend(x, node = node, method = method, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha, B = B, debug = debug,
    optimized = optimized)

}#LEARN.MB
