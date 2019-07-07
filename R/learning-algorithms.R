
# constraint-based learning algorithms.
bnlearn = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
    test = NULL, alpha = NULL, B = NULL, method = "gs", max.sx = NULL,
    debug = FALSE, optimized = FALSE, strict = FALSE, undirected = FALSE, 
    noise.levels = NULL) {

  reset.test.counter()

  res = NULL
  parallel = FALSE

  if (strict)
    warning("argument 'strict' is deprecated and will be removed in 2019.")
  if (optimized)
    warning("argument 'optimized' is deprecated and will be removed in 2019.")

  # check the data are there.
  data.info = check.data(x)
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
  # check size of the largest conditioning set in the independence tests.
  max.sx = check.largest.sx.set(max.sx, x)

  # check the cluster.
  if (!is.null(cluster)) {

    check.cluster(cluster)

    # enter in parallel mode.
    parallel = TRUE
    # set up the slave processes.
    slaves.setup(cluster)
    # disable debugging, the slaves do not cat() here.
    if (debug) {

      warning("disabling debugging output for parallel computing.")
      debug = FALSE

    }#THEN
    # disable backtracking, there's nothing to backtrack.
    optimized = FALSE

  }#THEN

  # sanitize whitelist and blacklist, if any.
  whitelist = build.whitelist(whitelist, nodes = names(x), data = x,
                algo = method, criterion = test)
  blacklist = build.blacklist(blacklist, whitelist, names(x), algo = method)
  # create the full blacklist incorporating model assumptions.
  full.blacklist = arcs.rbind(blacklist,
                     check.arcs.against.assumptions(NULL, x, test))
  full.blacklist = unique.arcs(full.blacklist, names(x))

  # call the right backend.
  if (method == "pc.stable") {

    local.structure =
      pc.stable.backend(x = x, whitelist = whitelist,
        blacklist = full.blacklist, test = test, alpha = alpha, B = B,
        max.sx = max.sx, debug = debug, cluster = cluster,
        complete = data.info$complete.nodes)

  }#THEN
  else if (method == "gs") {

    if (optimized) {

      local.structure =
        grow.shrink.optimized(x = x, whitelist = whitelist,
           blacklist = full.blacklist, test = test, alpha = alpha, B = B,
           max.sx = max.sx, strict = strict, debug = debug,
           complete = data.info$complete.nodes)

    }#THEN
    else {

      local.structure =
        grow.shrink(x = x, whitelist = whitelist, blacklist = full.blacklist,
          test = test, alpha = alpha, B = B, max.sx = max.sx,
          strict = strict, debug = debug, cluster = cluster,
          complete = data.info$complete.nodes)

    }#ELSE

  }#THEN
  else if (method == "iamb") {

    if (optimized) {

      local.structure =
        incremental.association.optimized(x = x, whitelist = whitelist,
          blacklist = full.blacklist, test = test, alpha = alpha, B = B,
          max.sx = max.sx, strict = strict, debug = debug,
          complete = data.info$complete.nodes)

    }#THEN
    else {

      local.structure =
        incremental.association(x = x, whitelist = whitelist,
          blacklist = full.blacklist, test = test, alpha = alpha, B = B,
          max.sx = max.sx, strict = strict, debug = debug,
          cluster = cluster, complete = data.info$complete.nodes)

    }#ELSE

  }#THEN
  else if (method == "fast.iamb") {

    if (optimized) {

      local.structure =
        fast.incremental.association.optimized(x = x, whitelist = whitelist,
          blacklist = full.blacklist, test = test, alpha = alpha, B = B,
          max.sx = max.sx, strict = strict, debug = debug,
          complete = data.info$complete.nodes)

    }#THEN
    else {

      local.structure =
        fast.incremental.association(x = x, whitelist = whitelist,
          blacklist = full.blacklist, test = test, alpha = alpha, B = B,
          max.sx = max.sx, strict = strict, debug = debug,
          cluster = cluster, complete = data.info$complete.nodes)

    }#ELSE

  }#THEN
  else if (method == "inter.iamb") {

    if (optimized) {

      local.structure =
        inter.incremental.association.optimized(x = x, whitelist = whitelist,
          blacklist = full.blacklist, test = test, alpha = alpha, B = B,
          max.sx = max.sx, strict = strict, debug = debug,
          complete = data.info$complete.nodes)

    }#THEN
    else {

      local.structure =
        inter.incremental.association(x = x, whitelist = whitelist,
          blacklist = full.blacklist, test = test, alpha = alpha, B = B,
          max.sx = max.sx, strict = strict, debug = debug,
          cluster = cluster, complete = data.info$complete.nodes)

    }#ELSE

  }#THEN
  else if (method == "mmpc") {

    if (optimized) {

      local.structure =
        maxmin.pc.optimized(x = x, whitelist = whitelist,
          blacklist = full.blacklist, test = test, alpha = alpha, B = B,
          max.sx = max.sx, strict = strict, debug = debug,
          complete = data.info$complete.nodes,
          noise.levels = noise.levels)

    }#THEN
    else {

      local.structure =
        maxmin.pc(x = x, whitelist = whitelist, blacklist = full.blacklist,
          test = test, alpha = alpha, B = B, strict = strict, debug = debug,
          max.sx = max.sx, cluster = cluster, complete = data.info$complete.nodes,
          noise.levels = noise.levels)

    }#ELSE

  }#THEN
  else if (method == "si.hiton.pc") {

    if (optimized) {

      local.structure =
        si.hiton.pc.optimized(x = x, whitelist = whitelist,
          blacklist = full.blacklist, test = test, alpha = alpha, B = B,
          max.sx = max.sx, strict = strict, debug = debug,
          complete = data.info$complete.nodes)

    }#THEN
    else {

      local.structure =
        si.hiton.pc.backend(x = x, whitelist = whitelist,
          blacklist = full.blacklist, test = test, alpha = alpha, B = B,
          max.sx = max.sx, strict = strict, debug = debug,
          cluster = cluster, complete = data.info$complete.nodes)

    }#ELSE

  }#THEN

  if (undirected) {

    # save the status of the learning algorithm.
    arcs = nbr2arcs(local.structure)
    learning = list(whitelist = whitelist, blacklist = blacklist,
      test = test, args = list(alpha = alpha), optimized = optimized,
      ntests = test.counter())

    # include also the number of permutations/bootstrap samples
    # if it makes sense.
    if (!is.null(B))
      learning$args$B = B

    res = list(learning = learning,
      nodes = cache.structure(names(x), arcs = arcs), arcs = arcs)

  }#THEN
  else {

    # recover some arc directions.
    res = second.principle(x = x, local.structure = local.structure,
            whitelist = whitelist, blacklist = full.blacklist, test = test,
            alpha = alpha, B = B, max.sx = max.sx, strict = strict,
            complete = data.info$complete.nodes, debug = debug)
    # return the user-specified blacklist, not the full one, as in other
    # classes of learning algorithms.
    res$learning["blacklist"] = list(blacklist)

  }#ELSE

  # add tests performed by the slaves to the test counter.
  if (parallel)
    res$learning$ntests = res$learning$ntests +
      sum(unlist(parallel::clusterEvalQ(cluster, test.counter())))
  # save the learning method used.
  res$learning$algo = method
  # save the 'optimized' flag.
  res$learning$optimized = optimized
  # save the 'undirected' flag.
  res$learning$undirected = undirected
  # save the maximum size of the tests' conditioning sets.
  res$learning$max.sx = max.sx
  # save arcs that are illegal according to parametric assumptions.
  res$learning$illegal = check.arcs.against.assumptions(NULL, x, test)

  invisible(structure(res, class = "bn"))

}#BNLEARN

# score-based learning algorithms.
greedy.search = function(x, start = NULL, whitelist = NULL, blacklist = NULL,
    score = NULL, heuristic = "hc", ..., optimized = TRUE, debug = FALSE) {

  # check the data are there.
  check.data(x)
  # check the algorithm.
  check.learning.algorithm(heuristic, class = "score")
  # check the score label.
  score = check.score(score, x)
  # check debug.
  check.logical(debug)

  # check unused arguments in misc.args.
  extra = list(...)
  misc.args = extra[names(extra) %in% method.extra.args[[heuristic]]]
  extra.args = extra[names(extra) %in% score.extra.args[[score]]]
  check.unused.args(extra, c(method.extra.args[[heuristic]], score.extra.args[[score]]))

  # expand and check the max.iter parameter (common to all algorithms).
  max.iter = check.max.iter(misc.args$max.iter)
  # expand and check the maxp parameter (common to all algorithms).
  maxp = check.maxp(misc.args$maxp, data = x)

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
  whitelist = build.whitelist(whitelist, nodes = names(x), data = x,
                algo = heuristic, criterion = score)
  blacklist = build.blacklist(blacklist, whitelist, names(x), algo = heuristic)
  # if there is no preseeded network, use an empty one.
  if (is.null(start))
    start = empty.graph(nodes = names(x))
  else {

    # check start's class.
    check.bn(start)
    # check the preseeded network against the data set.
    check.bn.vs.data(start, x)
    # check the preseeded network against the maximum number of parents.
    nparents = sapply(start$nodes, function(x) length(x$parents))
    if (any(nparents > maxp))
      stop("nodes ", paste(names(which(nparents > maxp)), collapse = " "),
        " have more than 'maxp' parents.")
    # check the preseeded network is valid for the model assumptions.
    check.arcs.against.assumptions(start$arcs, x, score)

  }#ELSE

  # apply the whitelist to the preseeded network; undirected arcs are allowed
  # but applied as directed to interoperate with bnlearn() in hybrid.search().
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
                 data = x, extra.args = extra.args, learning = TRUE)

  # reset the test counter.
  reset.test.counter()

  # call the right backend.
  if (heuristic == "hc") {

    res = hill.climbing(x = x, start = start, whitelist = whitelist,
      blacklist = blacklist, score = score, extra.args = extra.args,
      restart = restart, perturb = perturb, max.iter = max.iter,
      maxp = maxp, optimized = optimized, debug = debug)

  }#THEN
  else if (heuristic == "tabu"){

    res = tabu.search(x = x, start = start, whitelist = whitelist,
      blacklist = blacklist, score = score, extra.args = extra.args,
      max.iter = max.iter, optimized = optimized, tabu = tabu,
      maxp = maxp, debug = debug)

  }#THEN

  # set the metadata of the network in one stroke.
  res$learning = list(whitelist = whitelist, blacklist = blacklist,
    test = score, ntests = test.counter(),
    algo = heuristic, args = extra.args, optimized = optimized,
    illegal = check.arcs.against.assumptions(NULL, x, score))

  invisible(res)

}#GREEDY.SEARCH

# hybrid learning algorithms.
hybrid.search = function(x, whitelist = NULL, blacklist = NULL,
    restrict = "mmpc", maximize = "hc", restrict.args = list(),
    maximize.args = list(), debug = FALSE, noise.levels = NULL) {

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

  ## restrict phase.
  if (restrict %in% constraint.based.algorithms) {

    # merge the user-provided arguments with the defaults, making sure not to
    # overwrite critical arguments.
    critical.arguments = c("x", "method", "whitelist", "blacklist", "debug", "undirected", "noise.levels")
    named.arguments = names(formals(bnlearn))
    named.arguments = setdiff(named.arguments, critical.arguments)
    other.arguments = setdiff(names(restrict.args), named.arguments)
    check.unused.args(other.arguments, character(0))

    restrict.args[critical.arguments] =
      list(x, method = restrict, whitelist = whitelist, blacklist = blacklist,
           debug = debug, undirected = TRUE, noise.levels = noise.levels)

    rst = do.call("bnlearn", restrict.args)

  }#THEN
  else if (restrict %in% mim.based.algorithms) {

    check.unused.args(restrict.args, "mi")

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

  ## maximize phase.
  # merge the user-provided arguments with the defaults, making sure not to
  # overwrite critical arguments.
  critical.arguments = c("x", "start", "heuristic", "whitelist", "blacklist", "debug")
  named.arguments = names(formals(greedy.search))
  named.arguments = setdiff(named.arguments, critical.arguments)
  check.unused.args(intersect(critical.arguments, names(maximize.args)), character(0))

  maximize.args[critical.arguments] =
    list(x, start = NULL, heuristic = maximize, whitelist = whitelist, blacklist = constraints,
         debug = debug)

  res = do.call("greedy.search", maximize.args)

  # set the metadata of the network in one stroke.
  res$learning = list(whitelist = rst$learning$whitelist,
    blacklist = rst$learning$blacklist, test = res$learning$test,
    ntests = res$learning$ntests + rst$learning$ntests, algo = method,
    args = c(res$learning$args, rst$learning$args),
    optimized = res$learning$optimized || rst$learning$optimized,
    restrict = restrict, rstest = rst$learning$test, maximize = maximize,
    maxscore = res$learning$test,
    illegal = check.arcs.against.assumptions(NULL, x, rst$learning$test))

  invisible(res)

}#HYBRID.SEARCH

# learning algorithm based on the mutual information matrix.
mi.matrix = function(x, whitelist = NULL, blacklist = NULL, method, mi = NULL,
    debug = FALSE) {

  # check the data are there.
  check.data(x, allowed.types = c(discrete.data.types, continuous.data.types))
  # check the algorithm.
  check.learning.algorithm(method, class = "mim")
  # check debug.
  check.logical(debug)
  # check the label of the mutual information estimator.
  estimator = check.mi.estimator(mi, x)
  # sanitize whitelist and blacklist, if any.
  nodes = names(x)
  nnodes = length(nodes)

  whitelist = build.whitelist(whitelist, nodes = nodes, data = x,
                algo = method, criterion = estimator)
  blacklist = build.blacklist(blacklist, whitelist, nodes, algo = method)

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
    ntests = nnodes * (nnodes - 1) / 2, algo = method, undirected = TRUE,
    args = list(estimator = estimator))

  invisible(res)

}#MI.MATRIX

# learn the markov blanket of a single node.
mb.backend = function(x, target, method, whitelist = NULL, blacklist = NULL,
    start = NULL, test = NULL, alpha = 0.05, B = NULL, max.sx = NULL,
    debug = FALSE) {

  reset.test.counter()

  # check the data are there.
  data.info = check.data(x)
  # cache the node labels.
  nodes = names(x)
  # a valid node is needed.
  check.nodes(nodes = target, graph = nodes, max.nodes = 1)
  # check the algorithm.
  check.learning.algorithm(method, class = "markov.blanket")
  # check test labels.
  test = check.test(test, x)
  # check debug.
  check.logical(debug)
  # check alpha.
  alpha = check.alpha(alpha)
  # check B (the number of bootstrap/permutation samples).
  B = check.B(B, test)
  # check size of the largest conditioning set in the independence tests.
  max.sx = check.largest.sx.set(max.sx, x)

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

    blacklist = blacklist[blacklist %!in% whitelist]

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

    nodes = nodes[nodes %!in% blacklist]
    x = minimal.data.frame.column(x, nodes, drop = FALSE)
    x = minimal.data.frame(x)
    names(x) = nodes

  }#THEN

  # if all nodes are blacklisted, the markov blanket is obviously empty.
  if (all(nodes %in% c(target, blacklist)))
    return(character(0))

  # call the right backend.
  if (method == "gs") {

    mb = gs.markov.blanket(x = target, data = x, nodes = nodes, alpha = alpha,
           B = B, whitelist = whitelist, blacklist = NULL, start = start,
           backtracking = NULL, test = test, max.sx = max.sx,
           complete = data.info$complete.nodes, debug = debug)

  }#THEN
  else if (method == "iamb") {

    mb = ia.markov.blanket(x = target, data = x, nodes = nodes, alpha = alpha,
           B = B, whitelist = whitelist, blacklist = NULL, start = start,
           backtracking = NULL, test = test, max.sx = max.sx,
           complete = data.info$complete.nodes, debug = debug)

  }#THEN
  else if (method == "fast.iamb") {

    mb = fast.ia.markov.blanket(x = target, data = x, nodes = nodes,
           alpha = alpha, B = B, whitelist = whitelist, blacklist = NULL,
           start = start, backtracking = NULL, test = test,
		   complete = data.info$complete.nodes, max.sx = max.sx, debug = debug)

  }#THEN
  else if (method == "inter.iamb") {

    mb = inter.ia.markov.blanket(x = target, data = x, nodes = nodes,
           alpha = alpha, B = B, whitelist = whitelist, blacklist = NULL,
           start = start, backtracking = NULL, test = test, max.sx = max.sx,
           complete = data.info$complete.nodes, debug = debug)

  }#THEN

  return(mb)

}#MB.BACKEND

# learn the neighbourhood of a single node.
nbr.backend = function(x, target, method, whitelist = NULL, blacklist = NULL,
    test = NULL, alpha = 0.05, B = NULL, max.sx = NULL, debug = FALSE) {

  reset.test.counter()

  # check the data are there.
  data.info = check.data(x)
  # cache the node labels.
  nodes = names(x)
  # a valid node is needed.
  check.nodes(nodes = target, graph = nodes, max.nodes = 1)
  # check the algorithm.
  check.learning.algorithm(method, class = "neighbours")
  # check test labels.
  test = check.test(test, x)
  # check debug.
  check.logical(debug)
  # check alpha.
  alpha = check.alpha(alpha)
  # check B (the number of bootstrap/permutation samples).
  B = check.B(B, test)
  # check size of the largest conditioning set in the independence tests.
  max.sx = check.largest.sx.set(max.sx, x)

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

    blacklist = blacklist[blacklist %!in% whitelist]

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

    nodes = nodes[nodes %!in% blacklist]
    x = minimal.data.frame.column(x, nodes, drop = FALSE)
    x = minimal.data.frame(x)
    names(x) = nodes

  }#THEN

  # if all nodes are blacklisted, the neighbourhood is obviously empty.
  if (all(nodes %in% c(target, blacklist)))
    return(character(0))

  if (method %in% c("mmpc", "si.hiton.pc")) {

    # call the right backend, forward phase.
    if (method == "mmpc") {

      nbr = maxmin.pc.forward.phase(target, data = x, nodes = nodes,
             alpha = alpha, B = B, whitelist = whitelist, blacklist = blacklist,
             test = test, max.sx = max.sx, optimized = FALSE, debug = debug,
             backtracking = NULL, complete = data.info$complete.nodes)

    }#THEN
    else if (method == "si.hiton.pc") {

      nbr = si.hiton.pc.heuristic(target, data = x, nodes = nodes, alpha = alpha,
              B = B, whitelist = whitelist, blacklist = blacklist, test = test,
              max.sx = max.sx, backtracking = NULL, debug = debug,
              complete = data.info$complete.nodes)

    }#ELSE

    # this is the backward phase.
    nbr = neighbour(target, mb = structure(list(nbr), names = target), data = x,
            alpha = alpha, B = B, whitelist = whitelist, blacklist = blacklist,
            test = test, max.sx = max.sx, markov = FALSE, debug = debug,
            complete = data.info$complete.nodes)

    parents.and.children = nbr[["nbr"]]

  }#THEN
  else if (method == "pc.stable") {

    nbr = pc.stable.backend(x = x, whitelist = whitelist, blacklist = blacklist,
            test = test, alpha = alpha, B = B, max.sx = max.sx, debug = debug,
            complete = data.info$complete.nodes)

    parents.and.children = nbr[[target]]$nbr

  }#THEN

  return(parents.and.children)

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
  check.data(data, allowed.types = discrete.data.types)

  # check the explantory variables.
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
      check.nodes(explanatory, graph = explanatory)

  }#ELSE

  # check that at least one explanatory variable is provided.
  if (length(explanatory) == 0)
    stop("at least one explanatory variable is required.")
  # check that the training node is not included among the explanatory variables.
  if (training %in% explanatory)
    stop("node ", training, " is included in the model both as a training ",
         "and an explanatory variable.")
  # cache the whole node set.
  nodes = c(training, explanatory)
  # sanitize whitelist and blacklist, if any.
  if (method != "naive.bayes") {

    whitelist = build.whitelist(whitelist, nodes = nodes, data = data,
                  algo = method, criterion = "mi")
    blacklist = build.blacklist(blacklist, whitelist, nodes, algo = method)

    if (method == "tree.bayes") {

      # arcs to and from the training node cannot be whitelisted or blacklisted.
      if ((training %in% whitelist) || (training %in% blacklist))
        stop("blacklisting arcs to and from the training node is not allowed.")

    }#THEN

  }#THEN

  # sanitize method-specific arguments.
  extra.args = check.classifier.args(method = method, data = data, extra.args = expand,
                 training = training, explanatory = explanatory)

  if (method == "naive.bayes") {

    # naive bayes requires no test.
    ntests = 0
    # not test statistiic involved.
    test = "none"

    res = naive.bayes.backend(data = data, training = training,
            explanatory = explanatory)

  }#THEN
  else if (method == "tree.bayes") {

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
    test = test, ntests = ntests, algo = method, args = extra.args)
  # set the trainign variable, for use by predict() & co.
  res$learning$args$training = training

  invisible(res)

}#BAYESIAN.CLASSIFIER

