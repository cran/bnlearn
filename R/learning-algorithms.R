
# constraint-based learning algorithms.
bnlearn = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
    test = "mi", alpha = 0.05, B = NULL, method = "gs", debug = FALSE,
    optimized = TRUE, strict = TRUE, undirected = FALSE) {

  assign(".test.counter", 0, envir = .GlobalEnv)

  res = NULL
  cluster.aware = FALSE

  # check the data are there.
  check.data(x)
  # check the algorithm.
  check.learning.algorithm(method, class = "constraint")
  # check test labels.
  test = check.test(test, x)
  # check the logical flags (debug, strict, optimized, undirected).
  check.logical(debug)
  check.logical(strict)
  check.logical(optimized)
  check.logical(undirected)
  # check alpha.
  alpha = check.alpha(alpha)
  # check B (the number of bootstrap/permutation samples).
  B = check.B(B, test)

  # check the cluster.
  if (!is.null(cluster)) {

    check.cluster(cluster)

    # enter in cluster-aware mode.
    cluster.aware = TRUE
    # set the test counter in all the cluster nodes.
    clusterEvalQ(cluster, assign(".test.counter", 0, envir = .GlobalEnv))
    # disable debugging, the slaves do not cat() here.
    if (debug) {

      warning("disabling debugging output in cluster-aware mode.")
      debug = FALSE

    }#THEN

  }#THEN

  # sanitize whitelist and blacklist, if any.
  whitelist = build.whitelist(whitelist, names(x))
  blacklist = build.blacklist(blacklist, whitelist, names(x))

  # call the right backend.
  if (method == "gs") {

    if (cluster.aware) {

      mb = grow.shrink.cluster(x = x, cluster = cluster,
        whitelist = whitelist, blacklist = blacklist, test = test,
        alpha = alpha, B = B, strict = strict, debug = debug)

    }#THEN
    else if (optimized) {

      mb = grow.shrink.optimized(x = x, whitelist = whitelist,
        blacklist = blacklist, test = test, alpha = alpha, B = B,
        strict = strict, debug = debug)

    }#THEN
    else {

      mb = grow.shrink(x = x, whitelist = whitelist, blacklist = blacklist,
        test = test, alpha = alpha, B = B, strict = strict, debug = debug)

    }#ELSE

  }#THEN
  else if (method == "iamb") {

    if (cluster.aware) {

      mb = incremental.association.cluster(x = x, cluster = cluster,
        whitelist = whitelist, blacklist = blacklist, test = test,
        alpha = alpha, B = B, strict = strict, debug = debug)

    }#THEN
    else if (optimized) {

      mb = incremental.association.optimized(x = x, whitelist = whitelist,
        blacklist = blacklist, test = test, alpha = alpha, B = B,
        strict = strict, debug = debug)

    }#THEN
    else {

      mb = incremental.association(x = x, whitelist = whitelist,
        blacklist = blacklist, test = test, alpha = alpha, B = B,
        strict = strict, debug = debug)

    }#ELSE

  }#THEN
  else if (method == "fast.iamb") {

    if (cluster.aware) {

      mb = fast.incremental.association.cluster(x = x, cluster = cluster,
        whitelist = whitelist, blacklist = blacklist, test = test,
        alpha = alpha, B = B, strict = strict, debug = debug)

    }#THEN
    else if (optimized) {

      mb = fast.incremental.association.optimized(x = x, whitelist = whitelist,
        blacklist = blacklist, test = test, alpha = alpha, B = B,
        strict = strict, debug = debug)

    }#THEN
    else {

      mb = fast.incremental.association(x = x, whitelist = whitelist,
        blacklist = blacklist, test = test, alpha = alpha, B = B,
        strict = strict, debug = debug)

    }#ELSE

  }#THEN
  else if (method == "inter.iamb") {

    if (cluster.aware) {

      mb = inter.incremental.association.cluster(x = x, cluster = cluster,
        whitelist = whitelist, blacklist = blacklist, test = test,
        alpha = alpha, B = B, strict = strict, debug = debug)

    }#THEN
    else if (optimized) {

      mb = inter.incremental.association.optimized(x = x, whitelist = whitelist,
        blacklist = blacklist, test = test, alpha = alpha, B = B,
        strict = strict, debug = debug)

    }#THEN
    else {

      mb = inter.incremental.association(x = x, whitelist = whitelist,
        blacklist = blacklist, test = test, alpha = alpha, B = B,
        strict = strict, debug = debug)

    }#ELSE

  }#THEN
  else if (method == "mmpc") {

    if (cluster.aware) {

      mb = maxmin.pc.cluster(x = x, cluster = cluster, whitelist = whitelist,
        blacklist = blacklist, test = test, alpha = alpha, B = B,
        strict = strict, debug = debug)

    }#THEN
    else if (optimized) {

      mb = maxmin.pc.optimized(x = x, whitelist = whitelist, blacklist = blacklist,
        test = test, alpha = alpha, B = B, strict = strict, debug = debug)

    }#THEN
    else {

      mb = maxmin.pc(x = x, whitelist = whitelist, blacklist = blacklist,
        test = test, alpha = alpha, B = B, strict = strict, debug = debug)

    }#ELSE

  }#THEN
  else if (method == "si.hiton.pc") {

    if (cluster.aware) {

      mb = si.hiton.pc.cluster(x = x, cluster = cluster, whitelist = whitelist,
        blacklist = blacklist, test = test, alpha = alpha, B = B,
        strict = strict, debug = debug)

    }#THEN
    else if (optimized) {

      mb = si.hiton.pc.optimized(x = x, whitelist = whitelist, blacklist = blacklist,
        test = test, alpha = alpha, B = B, strict = strict, debug = debug)

    }#THEN
    else {

      mb = si.hiton.pc.backend(x = x, whitelist = whitelist, blacklist = blacklist,
        test = test, alpha = alpha, B = B, strict = strict, debug = debug)

    }#ELSE

  }#THEN

  if (undirected) {

    # save the status of the learning algorithm.
    arcs = nbr2arcs(mb)
    learning = list(whitelist = whitelist, blacklist = blacklist,
      test = test, args = list(alpha = alpha), optimized = optimized,
      ntests = get(".test.counter", envir = .GlobalEnv))

    # include also the number of permutations/bootstrap samples
    # if it makes sense.
    if (!is.null(B))
      learning$args$B = B

    res = list(learning = learning,
      nodes = cache.structure(names(mb), arcs = arcs), arcs = arcs)

  }#THEN
  else {

    # recover some arc directions.
    res = second.principle(x = x, mb = mb, whitelist = whitelist,
            blacklist = blacklist, test = test, alpha = alpha, B = B,
            strict = strict, debug = debug)

  }#ELSE

  # add tests performed by the slaves to the test counter.
  if (cluster.aware)
    res$learning$ntests = res$learning$ntests +
      sum(unlist(clusterEvalQ(cluster, get(".test.counter", envir = .GlobalEnv))))
  # save the learning method used.
  res$learning$algo = method
  # save the 'optimized' flag.
  res$learning$optimized = optimized

  invisible(structure(res, class = "bn"))

}#BNLEARN

# score-based learning algorithms.
greedy.search = function(x, start = NULL, whitelist = NULL, blacklist = NULL,
    score = "aic", heuristic = "hc", expand, optimized = FALSE, debug = FALSE) {

  # check the data are there.
  check.data(x)
  # check the algorithm.
  check.learning.algorithm(heuristic, class = "score")
  # check the score label.
  score = check.score(score, x)
  # check debug.
  check.logical(debug)

  # check unused arguments in misc.args.
  misc.args = expand[names(expand) %in% method.extra.args[[heuristic]]]
  extra.args = expand[names(expand) %in% score.extra.args[[score]]]
  check.unused.args(expand, c(method.extra.args[[heuristic]], score.extra.args[[score]]))

  # expand and check the max.iter parameter (common to all algorithm)
  max.iter = check.max.iter(misc.args$max.iter)

  if (heuristic == "hc") {

    # expand and check the number of random restarts.
    restart = check.restart(misc.args$restart)
    # expand and check the magnitude of the perturbation when random restarts
    # are effectively used.
    perturb = ifelse((restart > 0), check.perturb(misc.args$perturb), 0)

  }#THEN
  else if (heuristic == "tabu") {

    # expand and check the arguments related to the tabu list.
    tabu = check.tabu(misc.args$tabu)
    max.tabu = check.max.tabu(misc.args$max.tabu, tabu)

  }#THEN

  # sanitize whitelist and blacklist, if any.
  whitelist = build.whitelist(whitelist, names(x))
  blacklist = build.blacklist(blacklist, whitelist, names(x))
  # if there is no preseeded network, use an empty one.
  if (is.null(start))
    start = empty.graph(nodes = names(x))
  else {

    # check start's class.
    check.bn(start)
    # check the preseeded network against the data set.
    check.bn.vs.data(start, x)

  }#ELSE

  # apply the whitelist to the preseeded network.
  if (!is.null(whitelist)) {

    for (i in 1:nrow(whitelist))
      start$arcs = set.arc.direction(whitelist[i, "from"],
                       whitelist[i, "to"], start$arcs)

  }#THEN

  # apply the blacklist to the preseeded network.
  if (!is.null(blacklist)) {

    blacklisted = apply(start$arcs, 1, function(x){ is.blacklisted(blacklist, x) })
    start$arcs = start$arcs[!blacklisted, , drop = FALSE]

  }#THEN

  # be sure the graph structure is up to date.
  start$nodes = cache.structure(names(start$nodes), arcs = start$arcs)
  # no party if the graph is partially directed.
  if (is.pdag(start$arcs, names(start$nodes)))
    stop("the graph is only partially directed.")
  # check whether the graph is acyclic.
  if (!is.acyclic(arcs = start$arcs, nodes = names(start$nodes)))
    stop("the preseeded graph contains cycles.")

  # sanitize score-specific arguments.
  extra.args = check.score.args(score = score, network = start,
                 data = x, extra.args = extra.args)

  # create the test counter in .GlobalEnv.
  assign(".test.counter", 0, envir = .GlobalEnv)

  # call the right backend.
  if (heuristic == "hc") {

    res = hill.climbing(x = x, start = start, whitelist = whitelist,
      blacklist = blacklist, score = score, extra.args = extra.args,
      restart = restart, perturb = perturb, max.iter = max.iter,
      optimized = optimized, debug = debug)

  }#THEN
  else if (heuristic == "tabu"){

    res = tabu.search(x = x, start = start, whitelist = whitelist,
      blacklist = blacklist, score = score, extra.args = extra.args,
      max.iter = max.iter, optimized = optimized, tabu = tabu,
      debug = debug)

  }#THEN

  # set the metadata of the network in one stroke.
  res$learning = list(whitelist = whitelist, blacklist = blacklist,
    test = score, ntests = get(".test.counter", envir = .GlobalEnv),
    algo = heuristic, args = extra.args, optimized = optimized)

  invisible(res)

}#GREEDY.SEARCH

# hybrid learning algorithms.
hybrid.search = function(x, whitelist = NULL, blacklist = NULL,
    restrict = "mmpc", maximize = "hc", restrict.args = list(), score = NULL,
    maximize.args = list(), optimized = TRUE, debug = FALSE) {

  nodes = names(x)

  # check the restrict and maximize parameters.
  check.learning.algorithm(restrict, class = c("constraint", "mim"))
  check.learning.algorithm(maximize, class = "score")
  # choose the right method for the job.
  method = ifelse((restrict == "mmpc") && (maximize == "hc"), "mmhc", "rsmax2")

  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* restrict phase, using the", method.labels[restrict] ,"algorithm.\n")

  }#THEN

  # restrict phase
  if (restrict %in% constraint.based.algorithms) {

    rst = bnlearn(x, cluster = NULL, whitelist = whitelist, blacklist = blacklist,
            test = restrict.args$test, alpha = restrict.args$alpha,
            B = restrict.args$B, method = restrict, debug = debug,
            optimized = optimized, strict = restrict.args$strict, undirected = TRUE)

  }#THEN
  else if (restrict %in% mim.based.algorithms) {

    rst = mi.matrix(x, whitelist = whitelist, blacklist = blacklist,
            method = restrict, mi = restrict.args$mi, debug = debug)

  }#THEN

  # transform the constraints learned during the restrict phase in a blacklist
  # which will be used in the maximize phase.
  constraints = arcs.to.be.added(rst$arcs, nodes, whitelist = rst$learning$blacklist)

  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* maximize phase, using the", method.labels[maximize] ,"algorithm.\n")

  }#THEN

  # maximize phase
  res = greedy.search(x, start = NULL, whitelist = whitelist, blacklist = constraints,
          score = score, heuristic = maximize, expand = maximize.args,
          optimized = optimized, debug = debug)

  # set the metadata of the network in one stroke.
  res$learning = list(whitelist = rst$learning$whitelist,
    blacklist = rst$learning$blacklist, test = res$learning$test,
    ntests = res$learning$ntests + rst$learning$ntests, algo = method,
    args = c(res$learning$args, rst$learning$args), optimized = optimized,
    restrict = restrict, rstest = rst$learning$test, maximize = maximize,
    maxscore = res$learning$test)

  invisible(res)

}#HYBRID.SEARCH

# learning algorithm based on the mutual information matrix.
mi.matrix = function(x, whitelist = NULL, blacklist = NULL, method, mi = NULL,
    debug = FALSE) {

  # check the data are there.
  check.data(x)
  # check the algorithm.
  check.learning.algorithm(method, class = "mim")
  # check debug.
  check.logical(debug)
  # check the label of the mutual information estimator.
  estimator = check.mi.estimator(mi, x)
  # sanitize whitelist and blacklist, if any.
  nodes = names(x)
  nnodes = length(nodes)
  whitelist = build.whitelist(whitelist, nodes)
  blacklist = build.blacklist(blacklist, whitelist, nodes)

  if (method == "aracne") {

    arcs = aracne.backend(x = x, estimator = match(estimator, available.mi),
             whitelist = whitelist, blacklist = blacklist, debug = debug)

  }#THEN
  else if (method == "chow.liu") {

    # check whether any node has all incident arcs blacklisted; if so it's
    # simply not possible to learn a tree spanning all the nodes.
    culprit = names(which(table(blacklist) == 2 * (ncol(x) - 1)))

    if (length(culprit) > 0)
      stop("all arcs incident on nodes '", culprit, "' are blacklisted.")

    arcs = chow.liu.backend(x = x, nodes = nodes,
             estimator = match(estimator, available.mi), whitelist = whitelist,
             blacklist = blacklist, conditional = NULL, debug = debug)

  }#THEN

  res = empty.graph(nodes)
  # update the arcs of the network.
  res$arcs = arcs
  # update the network structure.
  res$nodes = cache.structure(nodes, arcs = arcs)
  # set the metadata of the network in one stroke.
  res$learning = list(whitelist = whitelist, blacklist = blacklist,
    test = as.character(mi.estimator.tests[estimator]), alpha = 0.05,
    ntests = nnodes * (nnodes - 1) / 2, algo = method,
    args = list(estimator = estimator), optimized = NULL)

  invisible(res)

}#MI.MATRIX

# learn the markov blanket of a single node.
mb.backend = function(x, target, method, whitelist = NULL, blacklist = NULL,
    start = NULL, test = NULL, alpha = 0.05, B = NULL, debug = FALSE,
    optimized = TRUE) {

  assign(".test.counter", 0, envir = .GlobalEnv)

  # check the data are there.
  check.data(x)
  # cache the node labels.
  nodes = names(x)
  # a valid node is needed.
  check.nodes(nodes = target, graph = nodes, max.nodes = 1)
  # check the algorithm.
  check.learning.algorithm(method, class = "markov.blanket")
  # check test labels.
  test = check.test(test, x)
  # check the logical flags (debug, optimized).
  check.logical(debug)
  check.logical(optimized)
  # check alpha.
  alpha = check.alpha(alpha)
  # check B (the number of bootstrap/permutation samples).
  B = check.B(B, test)

  # check the initial status of the markov blanket.
  if (!is.null(start)) {

    # must be made up of valid node labels.
    check.nodes(nodes = start, graph = nodes[nodes != target])

  }#THEN
  else {

    start = character(0)

  }#ELSE

  # sanitize and rework the whitelist.
  if (!is.null(whitelist)) {

    # target variable in the whitelist does not make much sense.
    if (target %in% whitelist)
      warning("target variable in the whitelist.")
    # check the labels of the whitelisted nodes.
    check.nodes(nodes = whitelist, graph = nodes[nodes != target])

    whitelist = matrix(c(rep(target, length(whitelist)), whitelist), ncol = 2)
    whitelist = arcs.rbind(whitelist, whitelist, reverse2 = TRUE)

  }#THEN

  # remove whitelisted nodes from the blacklist.
  if (!is.null(whitelist) && !is.null(blacklist)) {

    blacklist = blacklist[!(blacklist %in% whitelist)]

    if (length(blacklist) == 0)
      blacklist = NULL

  }#THEN

  # sanitize the blacklist, and drop the variables from the data.
  if (!is.null(blacklist)) {

    # target variable in the blacklist is plainly wrong.
    if (target %in% blacklist)
      stop("target variable in the blacklist.")
    # check the labels of the blacklisted nodes.
    check.nodes(nodes = blacklist, graph = nodes[nodes != target])

    nodes = nodes[!(nodes %in% blacklist)]
    x = minimal.data.frame.column(x, nodes, drop = FALSE)
    x = minimal.data.frame(x)
    names(x) = nodes

  }#THEN

  # call the right backend.
  if (method == "gs") {

    mb = gs.markov.blanket(x = target, data = x, nodes = nodes, alpha = alpha,
           B = B, whitelist = whitelist, blacklist = NULL, start = start,
           backtracking = NULL, test = test, debug = debug)

  }#THEN
  else if (method == "iamb") {

    mb = ia.markov.blanket(x = target, data = x, nodes = nodes, alpha = alpha,
           B = B, whitelist = whitelist, blacklist = NULL, start = start,
           backtracking = NULL, test = test, debug = debug)

  }#THEN
  else if (method == "fast.iamb") {

    mb = fast.ia.markov.blanket(x = target, data = x, nodes = nodes,
           alpha = alpha, B = B, whitelist = whitelist, blacklist = NULL,
           start = start, backtracking = NULL, test = test, debug = debug)

  }#THEN
  else if (method == "inter.iamb") {

    mb = inter.ia.markov.blanket(x = target, data = x, nodes = nodes, alpha = alpha,
           B = B, whitelist = whitelist, blacklist = NULL, start = start,
           backtracking = NULL, test = test, debug = debug)

  }#THEN

  return(mb)

}#MB.BACKEND

# learn the neighbourhood of a single node.
nbr.backend = function(x, target, method, whitelist = NULL, blacklist = NULL,
    test = NULL, alpha = 0.05, B = NULL, debug = FALSE, optimized = TRUE) {

  assign(".test.counter", 0, envir = .GlobalEnv)

  # check the data are there.
  check.data(x)
  # cache the node labels.
  nodes = names(x)
  # a valid node is needed.
  check.nodes(nodes = target, graph = nodes, max.nodes = 1)
  # check the algorithm.
  check.learning.algorithm(method, class = "neighbours")
  # check test labels.
  test = check.test(test, x)
  # check the logical flags (debug, optimized).
  check.logical(debug)
  check.logical(optimized)
  # check alpha.
  alpha = check.alpha(alpha)
  # check B (the number of bootstrap/permutation samples).
  B = check.B(B, test)

  # sanitize and rework the whitelist.
  if (!is.null(whitelist)) {

    # target variable in the whitelist does not make much sense.
    if (target %in% whitelist)
      warning("target variable in the whitelist.")
    # check the labels of the whitelisted nodes.
    check.nodes(nodes = whitelist, graph = nodes[nodes != target])

    whitelist = matrix(c(rep(target, length(whitelist)), whitelist), ncol = 2)
    whitelist = arcs.rbind(whitelist, whitelist, reverse2 = TRUE)

  }#THEN

  # remove whitelisted nodes from the blacklist.
  if (!is.null(whitelist) && !is.null(blacklist)) {

    blacklist = blacklist[!(blacklist %in% whitelist)]

    if (length(blacklist) == 0)
      blacklist = NULL

  }#THEN

  # sanitize and rework the blacklist.
  if (!is.null(blacklist)) {

    # target variable in the blacklist is plainly wrong.
    if (target %in% blacklist)
      stop("target variable in the blacklist.")
    # check the labels of the blacklisted nodes.
    check.nodes(nodes = blacklist, graph = nodes[nodes != target])

    nodes = nodes[!(nodes %in% blacklist)]
    x = minimal.data.frame.column(x, nodes, drop = FALSE)
    x = minimal.data.frame(x)
    names(x) = nodes

  }#THEN

  # call the right backend, forward phase.
  if (method == "mmpc") {

    nbr = maxmin.pc.forward.phase(target, data = x, nodes = nodes, 
           alpha = alpha, B = B, whitelist = whitelist, blacklist = blacklist,
           test = test, optimized = optimized, debug = debug)

  }#THEN
  else if (method == "si.hiton.pc") {

    nbr = si.hiton.pc.heuristic(target, data = x, nodes = nodes, alpha = alpha,
            B = B, whitelist = whitelist, blacklist = blacklist, test = test,
            optimized = optimized, debug = debug) 

  }#ELSE

  # this is the backward phase.
  nbr = neighbour(target, mb = structure(list(nbr), names = target), data = x, 
          alpha = alpha, B = B, whitelist = whitelist, blacklist = blacklist,
          test = test, markov = FALSE, debug = debug)

  return(nbr[["nbr"]])

}#NBR.BACKEND

# baeysian network classifiers.
bayesian.classifier = function(data, method, training, explanatory, whitelist,
    blacklist, expand, debug = FALSE) {

  # check debug.
  check.logical(debug)
  # check the learning algorithm.
  check.learning.algorithm(method, class = "classifier")
  # check the training node (the center of the star-shaped graph).
  check.nodes(training, max.nodes = 1)
  # check the data.
  if (method != "naive") {

      check.data(data)
    if (!is.data.discrete(data))
      stop("continuous data are not supported.")

  }#THEN

  # check the explantory variables (the points of the star-shaped graph).
  if (missing(data)) {

    check.nodes(explanatory)

  }#THEN
  else {

    vars = names(data)
    # check the label of the training variable.
    check.nodes(training, graph = vars, max.nodes = 1)
    # check the labels of the explanatory variables.
    if (missing(explanatory))
      explanatory = vars[vars != training]
    else
      check.nodes(explanatory, graph = nodes)

  }#ELSE

  # check that the training node is not included among the explanatory variables.
  if (training %in% explanatory)
    stop("node ", training, " is included in the model both as a training ",
         "and an explanatory variable.")
  # cache the whole node set.
  nodes = c(training, explanatory)
  # sanitize whitelist and blacklist, if any.
  if (method != "naive") {

    whitelist = build.whitelist(whitelist, explanatory)
    blacklist = build.blacklist(blacklist, whitelist, explanatory)

  }#THEN

  # sanitize method-specific arguments.
  extra.args = check.classifier.args(method = method, data = data, extra.args = expand,
                 training = training, explanatory = explanatory)

  if (method == "naive") {

    # naive bayes requires no test.
    ntests = 0
    # not test statistiic involved.
    test = "none"

    res = naive.bayes.backend(data = data, training = training,
            explanatory = explanatory)

  }#THEN
  else if (method == "tan") {

    # tan gets its tests from the chow-liu algorithm.
    ntests = length(explanatory) * (length(explanatory) - 1)/2
    # same for the test
    test = as.character(mi.estimator.tests[extra.args$estimator])

    res = tan.backend(data = data, training = training, explanatory = explanatory,
            whitelist = whitelist, blacklist = blacklist, mi = extra.args$estimator,
            root = extra.args$root, debug = debug)

  }#THEN

  # set the learning algorithm.
  res$learning$algo = method
  # set the metadata of the network in one stroke.
  res$learning = list(whitelist = whitelist, blacklist = blacklist,
    test = test, ntests = ntests, algo = method, args = extra.args,
    optimized = NULL)
  # set the trainign variable, for use by predict() & co.
  res$learning$args$training = training

  invisible(res)

}#BAYESIAN.CLASSIFIER

