
# constraint-based learning algorithms.
bnlearn = function(data, cluster = NULL, whitelist, blacklist, test = NULL,
    alpha = NULL, extra.args = list(), algorithm, max.sx = NULL, debug = FALSE,
    undirected = FALSE) {

  reset.test.counter()

  parallel = FALSE

  # check the data are there.
  data = check.data(data, allow.missing = TRUE)
  nodes = names(data)
  # check the constraint-based structure learning algorithm.
  check.learning.algorithm(algorithm, class = "constraint")
  # check test labels.
  test = check.test(test, data = data)
  # check the logical flags (debug, undirected).
  check.logical(debug)
  check.logical(undirected)
  # check alpha.
  alpha = check.alpha(alpha)
  # check the optional arguments to the test.
  extra.args = check.test.args(test = test, extra.args = extra.args, data = data)
  # check size of the largest conditioning set in the independence tests.
  max.sx = check.largest.sx.set(max.sx, data)

  # check the cluster.
  cluster = check.cluster(cluster)

  if (!is.null(cluster)) {

    # enter in parallel mode.
    parallel = TRUE
    # set up the slave processes.
    slaves.setup(cluster)
    # disable debugging, the slaves do not cat() here.
    if (debug) {

      warning("disabling debugging output for parallel computing.")
      debug = FALSE

    }#THEN

  }#THEN

  # sanitize whitelist and blacklist, if any.
  whitelist = build.whitelist(whitelist, nodes = nodes, data = data,
                algo = algorithm, criterion = test)
  blacklist = build.blacklist(blacklist, whitelist, nodes, algo = algorithm)
  # create the full blacklist incorporating model assumptions.
  full.blacklist = arcs.rbind(blacklist, list.illegal.arcs(nodes, data, test))
  full.blacklist = arcs.unique(full.blacklist, nodes)

  # call the right backend.
  if (algorithm == "pc.stable") {

    local.structure =
      pc.stable.backend(data = data, whitelist = whitelist,
        blacklist = full.blacklist, test = test, alpha = alpha,
        extra.args = extra.args, max.sx = max.sx, debug = debug,
        cluster = cluster)

  }#THEN
  else if (algorithm == "gs") {

    local.structure =
      grow.shrink(data = data, whitelist = whitelist,
        blacklist = full.blacklist, test = test, alpha = alpha,
        extra.args = extra.args, max.sx = max.sx, debug = debug,
        cluster = cluster)

  }#THEN
  else if (algorithm == "iamb") {

    local.structure =
      incremental.association(data = data, whitelist = whitelist,
        blacklist = full.blacklist, test = test, alpha = alpha,
        extra.args = extra.args, max.sx = max.sx, debug = debug,
        cluster = cluster)

  }#THEN
  else if (algorithm == "fast.iamb") {

    warning("fast.iamb() is deprecated and will be removed in 2027.")

    local.structure =
      fast.incremental.association(data = data, whitelist = whitelist,
        blacklist = full.blacklist, test = test, alpha = alpha,
        extra.args = extra.args, max.sx = max.sx, debug = debug,
        cluster = cluster)

  }#THEN
  else if (algorithm == "inter.iamb") {

    local.structure =
      inter.incremental.association(data = data, whitelist = whitelist,
        blacklist = full.blacklist, test = test, alpha = alpha,
        extra.args = extra.args, max.sx = max.sx, debug = debug,
        cluster = cluster)

  }#THEN
  else if (algorithm == "iamb.fdr") {

    local.structure =
      incremental.association.fdr(data = data, whitelist = whitelist,
        blacklist = full.blacklist, test = test, alpha = alpha,
        extra.args = extra.args, max.sx = max.sx, debug = debug,
        cluster = cluster)

  }#THEN
  else if (algorithm == "mmpc") {

    local.structure =
      maxmin.pc(data = data, whitelist = whitelist, blacklist = full.blacklist,
        test = test, alpha = alpha, extra.args = extra.args, debug = debug,
        max.sx = max.sx, cluster = cluster)

  }#THEN
  else if (algorithm == "si.hiton.pc") {

    local.structure =
      si.hiton.pc.backend(data = data, whitelist = whitelist,
        blacklist = full.blacklist, test = test, alpha = alpha,
        extra.args = extra.args, max.sx = max.sx, debug = debug,
        cluster = cluster)

  }#THEN
  else if (algorithm == "hpc") {

      local.structure =
        hybrid.pc.backend(data = data, whitelist = whitelist,
          blacklist = full.blacklist, test = test, alpha = alpha,
          extra.args = extra.args, max.sx = max.sx, debug = debug,
          cluster = cluster)

  }#THEN

  # recover some arc directions, or not.
  if (undirected) {

    arcs = sets2arcs(local.structure, names(local.structure))

  }#THEN
  else {

    arcs = learn.arc.directions(data = data, local.structure = local.structure,
             whitelist = whitelist, blacklist = full.blacklist, test = test,
             alpha = alpha, extra.args = extra.args, max.sx = max.sx,
             debug = debug)

  }#ELSE

  # group all the test arguments.
  all.args = c(list(alpha = alpha), extra.args)
  # join them with the relevant counters and learning arguments.
  learning = list(whitelist = whitelist, blacklist = blacklist, test = test,
                  args = all.args, ntests = test.counter(),
                  algo = algorithm, undirected = undirected, max.sx = max.sx)
  # add the list of illegal arcs, if we have one.
  learning$illegal = list.illegal.arcs(nodes, data, test)
  # add tests performed by the slaves to the test counter.
  if (parallel)
    learning$ntests = learning$ntests +
      sum(unlist(parallel::clusterEvalQ(cluster, test.counter())))

  pdag = structure(list(learning = learning,
                        nodes = cache.structure(nodes, arcs = arcs),
                        arcs = arcs),
                   class = "bn")

  invisible(pdag)

}#BNLEARN

# score-based learning algorithms.
greedy.search = function(data, start = NULL, whitelist, blacklist, score = NULL,
    algorithm, ..., optimized = TRUE, debug = FALSE) {

  # check the data are there.
  data = check.data(data, allow.missing = TRUE)
  nodes = names(data)
  # check the score-based structure learning algorithm.
  check.learning.algorithm(algorithm, class = "score")
  # check the score label.
  score = check.score(score, data = data)

  check.logical(debug)

  # check unused arguments in misc.args.
  extra = list(...)
  misc.args = extra[names(extra) %in% learning.extra.args[[algorithm]]]
  extra.args = extra[names(extra) %in% score.extra.args[[score]]]
  check.unused.args(extra, c(learning.extra.args[[algorithm]],
                             score.extra.args[[score]]))

  # expand and check the max.iter argument (common to all algorithms).
  max.iter = check.max.iter(misc.args$max.iter)
  # expand and check the maxp argument (common to all algorithms).
  maxp = check.maxp(misc.args$maxp, nnodes = ncol(data))

  if (algorithm == "hc") {

    # expand and check the number of random restarts.
    restart = check.restart(misc.args$restart)
    # expand and check the magnitude of the perturbation when random restarts
    # are effectively used.
    perturb = ifelse((restart > 0), check.perturb(misc.args$perturb), 0)

  }#THEN
  else if (algorithm == "tabu") {

    # expand and check the arguments related to the tabu list.
    tabu = check.tabu(misc.args$tabu)
    max.tabu = check.max.tabu(misc.args$max.tabu, tabu)

  }#THEN

  # sanitize whitelist and blacklist, if any.
  whitelist = build.whitelist(whitelist, nodes = nodes, data = data,
                algo = algorithm, criterion = score)
  blacklist = build.blacklist(blacklist, whitelist, nodes, algo = algorithm)
  # if there is no preseeded network, use an empty one.
  if (is.null(start))
    start = empty.graph(nodes = nodes)
  else {

    check.bn(start)
    # check the preseeded network against the data set.
    check.bn.vs.data(start, data)
    # check whether the network is valid for the score.
    check.arcs.against.assumptions(start$arcs, data, score)
    # check the preseeded network against the maximum number of parents.
    nparents = sapply(start$nodes, function(x) length(x$parents))
    if (any(nparents > maxp))
      stop("nodes ", paste(names(which(nparents > maxp)), collapse = " "),
        " have more than 'maxp' parents.")
    # check the preseeded network is valid for the model assumptions.
    check.arcs.against.assumptions(start$arcs, data, score)

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
  if (!is.completely.directed(start))
    stop("the graph is only partially directed.")
  # check whether the graph is acyclic.
  if (!is.acyclic(arcs = start$arcs, nodes = names(start$nodes)))
    stop("the preseeded graph contains cycles.")

  # sanitize score-specific arguments.
  extra.args = check.score.args(score = score, network = start,
                 data = data, extra.args = extra.args, learning = TRUE)

  # reset the test counter.
  reset.test.counter()

  # call the right backend.
  if (algorithm == "hc") {

    dag = hill.climbing(x = data, start = start, whitelist = whitelist,
      blacklist = blacklist, score = score, extra.args = extra.args,
      restart = restart, perturb = perturb, max.iter = max.iter,
      maxp = maxp, optimized = optimized, debug = debug)

  }#THEN
  else if (algorithm == "tabu"){

    dag = tabu.search(x = data, start = start, whitelist = whitelist,
      blacklist = blacklist, score = score, extra.args = extra.args,
      max.iter = max.iter, optimized = optimized, tabu = tabu,
      maxp = maxp, debug = debug)

  }#THEN

  # set the metadata of the network in one stroke.
  dag$learning = list(whitelist = whitelist, blacklist = blacklist,
    test = score, ntests = test.counter(),
    algo = algorithm, args = extra.args, optimized = optimized,
    illegal = list.illegal.arcs(nodes, data = data, criterion = score))

  invisible(dag)

}#GREEDY.SEARCH

# hybrid learning algorithms.
hybrid.search = function(data, whitelist, blacklist, restrict, maximize,
    restrict.args = list(), maximize.args = list(), debug = FALSE) {

  nodes = names(data)

  # check the restrict and maximize arguments.
  check.learning.algorithm(restrict, class = c("constraint", "mim"))
  check.learning.algorithm(maximize, class = "score")
  # choose the right algorithm for the job.
  if ((restrict == "mmpc") && (maximize == "hc"))
    algorithm = "mmhc"
  else if ((restrict == "hpc") && (maximize == "hc"))
    algorithm = "h2pc"
  else
   algorithm = "rsmax2"

  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* restrict phase, using the", learning.labels[restrict] ,"algorithm.\n")

  }#THEN

  # stub whitelists and blacklists before referencing them.
  if (missing(blacklist))
    blacklist = NULL
  if (missing(whitelist))
    whitelist = NULL

  ## restrict phase.
  if (restrict %in% constraint.based.algorithms) {

    # optional arguments to tests should be encapsulated into the extra.args
    # argument expected by bnlearn().
    if ("test" %in% names(restrict.args)) {

      test.label = restrict.args[["test"]]
      test.args = intersect(test.extra.args[[test.label]], names(restrict.args))
      extra.args = restrict.args[test.args]
      restrict.args[["extra.args"]] = extra.args
      restrict.args[test.args] = NULL

    }#THEN

    # merge the user-provided arguments with the defaults, making sure not to
    # overwrite critical arguments.
    critical.arguments =
      c("data", "algorithm", "whitelist", "blacklist", "debug", "undirected")
    named.arguments = names(formals(bnlearn))
    named.arguments = setdiff(named.arguments, critical.arguments)
    other.arguments = setdiff(names(restrict.args), named.arguments)
    check.unused.args(other.arguments, character(0))

    restrict.args[critical.arguments] =
      list(data, algorithm = restrict, whitelist = whitelist,
           blacklist = blacklist, debug = debug, undirected = TRUE)

    rst = do.call("bnlearn", restrict.args)

  }#THEN
  else if (restrict %in% mim.based.algorithms) {

    # warn about and remove unused arguments.
    restrict.args = check.unused.args(restrict.args, "mi")

    rst = mi.matrix(data, whitelist = whitelist, blacklist = blacklist,
            algorithm = restrict, mi = restrict.args$mi, debug = debug)

  }#THEN

  # transform the constraints learned during the restrict phase in a blacklist
  # which will be used in the maximize phase.
  constraints = arcs.to.be.added(rst$arcs, nodes, whitelist = rst$learning$blacklist)

  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* maximize phase, using the", learning.labels[maximize] ,"algorithm.\n")

  }#THEN

  ## maximize phase.
  # merge the user-provided arguments with the defaults, making sure not to
  # overwrite critical arguments.
  critical.arguments =
    c("data", "start", "algorithm", "whitelist", "blacklist", "debug")
  named.arguments = names(formals(greedy.search))
  named.arguments = setdiff(named.arguments, critical.arguments)
  check.unused.args(intersect(critical.arguments, names(maximize.args)),
                    character(0))

  # constraint-based algorithms allow whitelisting an arc in both directions,
  # but score-based algorithms do not; if the whitelist contains undirected arcs
  # take it as graph and use its consistent extension as a compromise.
  if (any(which.undirected(whitelist, nodes = nodes)))
    whitelist = .Call(call_pdag_extension,
                      arcs = whitelist,
                      nodes = nodes,
                      debug = FALSE)

  maximize.args[critical.arguments] =
    list(data, start = NULL, algorithm = maximize, whitelist = whitelist,
         blacklist = constraints, debug = debug)

  dag = do.call("greedy.search", maximize.args)

  # set the metadata of the network in one stroke.
  dag$learning = list(whitelist = rst$learning$whitelist,
    blacklist = rst$learning$blacklist, test = dag$learning$test,
    ntests = dag$learning$ntests + rst$learning$ntests, algo = algorithm,
    args = c(dag$learning$args, rst$learning$args),
    optimized = dag$learning$optimized,
    restrict = restrict, rstest = rst$learning$test, maximize = maximize,
    maxscore = dag$learning$test,
    illegal = list.illegal.arcs(nodes, data = data, criterion = dag$learning$test))

  invisible(dag)

}#HYBRID.SEARCH

# learning algorithm based on the mutual information matrix.
mi.matrix = function(data, whitelist, blacklist, algorithm, mi = NULL,
    debug = FALSE) {

  # check the data are there.
  data = check.data(data, allow.missing = TRUE, stop.if.all.missing = TRUE,
        allowed.types = c(discrete.data.types, continuous.data.types))
  nodes = names(data)
  nnodes = length(nodes)
  # check the structure learning algorithm.
  check.learning.algorithm(algorithm, class = "mim")
  # check the label of the mutual information estimator.
  estimator = check.mi.estimator(mi, data)
  # sanitize whitelist and blacklist, if any.
  whitelist = build.whitelist(whitelist, nodes = nodes, data = data,
                algo = algorithm, criterion = estimator)
  blacklist = build.blacklist(blacklist, whitelist, nodes, algo = algorithm)

  check.logical(debug)

  if (algorithm == "aracne") {

    arcs = aracne.backend(data = data, estimator = estimator,
             whitelist = whitelist, blacklist = blacklist, debug = debug)

  }#THEN
  else if (algorithm == "chow.liu") {

    # check whether any node has all incident arcs blacklisted; if so it's
    # simply not possible to learn a tree spanning all the nodes.
    culprit = names(which(table(blacklist) == 2 * (ncol(data) - 1)))

    if (length(culprit) > 0)
      stop("all arcs incident on nodes '", culprit, "' are blacklisted.")

    arcs = chow.liu.backend(data = data, nodes = nodes, estimator = estimator,
             whitelist = whitelist, blacklist = blacklist, conditional = NULL,
             debug = debug)

  }#THEN

  ugraph = empty.graph(nodes)
  # update the arcs of the network.
  ugraph$arcs = arcs
  # update the network structure.
  ugraph$nodes = cache.structure(nodes, arcs = arcs)
  # set the metadata of the network in one stroke.
  ugraph$learning = list(whitelist = whitelist, blacklist = blacklist,
    test = as.character(mi.estimator.tests[estimator]), alpha = 0.05,
    ntests = nnodes * (nnodes - 1) / 2, algo = algorithm, undirected = TRUE,
    args = list(estimator = estimator))

  invisible(ugraph)

}#MI.MATRIX

# learn the markov blanket of a single node.
mb.backend = function(data, target, algorithm, whitelist, blacklist,
    start = NULL, test = NULL, alpha = 0.05, extra.args = list(), max.sx = NULL,
    debug = FALSE) {

  reset.test.counter()

  # check the data are there.
  data = check.data(data, allow.missing = TRUE)
  nodes = names(data)
  # a valid node is needed.
  check.nodes(nodes = target, graph = nodes, max.nodes = 1)
  # check the Markov blanket detection algorithm.
  check.learning.algorithm(algorithm, class = "markov.blanket")
  # check test labels.
  test = check.test(test, data = data)

  check.logical(debug)
  # check alpha.
  alpha = check.alpha(alpha)
  # check the optional arguments to the test.
  extra.args = check.test.args(test = test, extra.args = extra.args, data = data)
  # check size of the largest conditioning set in the independence tests.
  max.sx = check.largest.sx.set(max.sx, data)

  # check the initial status of the markov blanket.
  if (!is.null(start)) {

    # must be made up of valid node labels.
    check.nodes(nodes = start, graph = nodes[nodes != target])

  }#THEN
  else {

    start = character(0)

  }#ELSE

  # stub whitelists and blacklists before referencing them.
  if (missing(blacklist))
    blacklist = NULL
  if (missing(whitelist))
    whitelist = NULL

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
    data = data[, nodes, drop = FALSE]

    # re-validate the data.
    data = check.data(data, allow.missing = TRUE)

  }#THEN

  # if all nodes are blacklisted, the markov blanket is obviously empty.
  if (all(nodes %in% c(target, blacklist)))
    return(character(0))

  # call the right backend.
  if (algorithm == "gs") {

    mb = gs.markov.blanket(target = target, data = data, nodes = nodes,
           alpha = alpha, extra.args = extra.args, whitelist = whitelist,
           blacklist = NULL, start = start, test = test, max.sx = max.sx,
           debug = debug)

  }#THEN
  else if (algorithm == "iamb") {

    mb = ia.markov.blanket(target = target, data = data, nodes = nodes,
           alpha = alpha, extra.args = extra.args, whitelist = whitelist,
           blacklist = NULL, start = start, test = test, max.sx = max.sx,
           debug = debug)

  }#THEN
  else if (algorithm == "fast.iamb") {

    warning("fast.iamb() is deprecated and will be removed in 2027.")

    mb = fast.ia.markov.blanket(target = target, data = data, nodes = nodes,
           alpha = alpha, extra.args = extra.args, whitelist = whitelist,
           blacklist = NULL, start = start, test = test, max.sx = max.sx,
           debug = debug)

  }#THEN
  else if (algorithm == "inter.iamb") {

    mb = inter.ia.markov.blanket(target = target, data = data, nodes = nodes,
           alpha = alpha, extra.args = extra.args, whitelist = whitelist,
           blacklist = NULL, start = start, test = test, max.sx = max.sx,
           debug = debug)

  }#THEN
  else if (algorithm == "iamb.fdr") {

    mb = ia.fdr.markov.blanket(target = target, data = data, nodes = nodes,
           alpha = alpha, extra.args = extra.args, whitelist = whitelist,
           blacklist = NULL, start = start, test = test, max.sx = max.sx,
           debug = debug)

  }#THEN

  return(mb)

}#MB.BACKEND

# learn the neighbourhood of a single node.
nbr.backend = function(data, target, algorithm, whitelist, blacklist,
    test = NULL, alpha = 0.05, extra.args = list(), max.sx = NULL,
    debug = FALSE) {

  reset.test.counter()

  # check the data are there.
  data = check.data(data, allow.missing = TRUE)
  nodes = names(data)
  # a valid node is needed.
  check.nodes(nodes = target, graph = nodes, max.nodes = 1)
  # check the feature selection algorithm.
  check.learning.algorithm(algorithm, class = "neighbours")
  # check test labels.
  test = check.test(test, data = data)

  check.logical(debug)
  # check alpha.
  alpha = check.alpha(alpha)
  # check the optional arguments to the test.
  extra.args = check.test.args(test = test, extra.args = extra.args, data = data)
  # check size of the largest conditioning set in the independence tests.
  max.sx = check.largest.sx.set(max.sx, data)

  # stub whitelists and blacklists before referencing them.
  if (missing(blacklist))
    blacklist = NULL
  if (missing(whitelist))
    whitelist = NULL

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
    data = data[, nodes, drop = FALSE]
    names(data) = nodes

  }#THEN

  # if all nodes are blacklisted, the neighbourhood is obviously empty.
  if (all(nodes %in% c(target, blacklist)))
    return(character(0))

  if (algorithm %in% c("mmpc", "si.hiton.pc")) {

    # call the right backend, forward phase.
    if (algorithm == "mmpc") {

      nbr = maxmin.pc.forward.phase(target = target, data = data, nodes = nodes,
             alpha = alpha, extra.args = extra.args, whitelist = whitelist,
             blacklist = blacklist, test = test, max.sx = max.sx,
             debug = debug)

    }#THEN
    else if (algorithm == "si.hiton.pc") {

      nbr = si.hiton.pc.heuristic(target = target, data = data, nodes = nodes,
              alpha = alpha, extra.args = extra.args, whitelist = whitelist,
              blacklist = blacklist, test = test, max.sx = max.sx,
              debug = debug)

    }#ELSE

    # this is the backward phase.
    nbr = neighbour(target = target, mb = structure(list(nbr), names = target),
            data = data, alpha = alpha, extra.args = extra.args,
            whitelist = whitelist, blacklist = blacklist, test = test,
            max.sx = max.sx, markov = FALSE, debug = debug)

    parents.and.children = nbr[["nbr"]]

  }#THEN
  else if (algorithm == "hpc") {

    nbr = hybrid.pc.heuristic(target, data = data, nodes = nodes, alpha = alpha,
            extra.args = extra.args, whitelist = whitelist,
            blacklist = blacklist, test = test, max.sx = max.sx, debug = debug)

    parents.and.children = nbr$nbr

  }#THEN
  else if (algorithm == "pc.stable") {

    nbr = pc.stable.backend(data = data, whitelist = whitelist,
            blacklist = blacklist, test = test, alpha = alpha,
            extra.args = extra.args, max.sx = max.sx, debug = debug)

    parents.and.children = nbr[[target]]$nbr

  }#THEN

  return(parents.and.children)

}#NBR.BACKEND

# bayesian network classifiers.
bayesian.classifier = function(data, algorithm, training, explanatory,
    whitelist, blacklist, expand, debug = FALSE) {

  check.logical(debug)
  # check the learning algorithm.
  check.learning.algorithm(algorithm, class = "classifier")
  # check the training node (the center of the star-shaped graph).
  check.nodes(training, max.nodes = 1)
  # check the data.
  data = check.data(data, allowed.types = discrete.data.types,
           allow.missing = TRUE, stop.if.all.missing = TRUE)

  # check the explanatory variables.
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
  if (algorithm == "tree.bayes") {

    whitelist = build.whitelist(whitelist, nodes = nodes, data = data,
                  algo = "chow.liu", criterion = "mi")
    blacklist = build.blacklist(blacklist, whitelist, nodes, algo = "chow.liu")

    # arcs to and from the training node cannot be whitelisted or blacklisted.
    if ((training %in% whitelist) || (training %in% blacklist))
      stop("blacklisting arcs to and from the training node is not allowed.")

  }#THEN

  # sanitize algorithm-specific arguments.
  extra.args = check.classifier.args(algorithm = algorithm, data = data,
                 extra.args = expand, training = training,
                 explanatory = explanatory)

  if (algorithm == "naive.bayes") {

    # naive bayes requires no test.
    ntests = 0
    # not test statistic involved.
    test = "none"

    dag = naive.bayes.backend(data = data, training = training,
            explanatory = explanatory)

  }#THEN
  else if (algorithm == "tree.bayes") {

    # tan gets its tests from the chow-liu algorithm.
    ntests = length(explanatory) * (length(explanatory) - 1)/2
    # same for the test
    test = as.character(mi.estimator.tests[extra.args$estimator])

    dag = tan.backend(data = data, training = training,
            explanatory = explanatory, whitelist = whitelist,
            blacklist = blacklist, mi = extra.args$estimator,
            root = extra.args$root, debug = debug)

  }#THEN

  # set the learning algorithm.
  dag$learning$algo = algorithm
  # set the metadata of the network in one stroke.
  dag$learning = list(whitelist = whitelist, blacklist = blacklist,
    test = test, ntests = ntests, algo = algorithm, args = extra.args)
  # set the training variable, for use by predict() & co.
  dag$learning$args$training = training

  invisible(dag)

}#BAYESIAN.CLASSIFIER

# causal discovery based on the LiNGAM principle.
lingam.learners = function(data, cluster, algorithm, whitelist, blacklist, mi,
    maximize = "alasso", maximize.args = list(), debug = FALSE) {

  reset.test.counter()

  parallel = FALSE

  # check the data are there.
  data = check.data(data, allowed.types = "continuous", allow.missing = FALSE)
  nodes = names(data)
  # check the causal discovery algorithm.
  algorithm = check.learning.algorithm(algorithm, class = "lingam")
  # check the mutual information label.
  if (missing(mi)) {

    mi = "pwling"

  }#THEN
  else {

    check.label(mi, choices = c("pwling", "gkernel"),
      labels = c("Pairwise Non-Gauss.", "Gauss. Kernel"),
      argname = "mutual information estimator")

  }#ELSE
  # sanitize whitelist and blacklist, if any.
  whitelist = build.whitelist(whitelist, nodes = nodes, data = data,
                algo = algorithm, criterion = mi)
  blacklist = build.blacklist(blacklist, whitelist, nodes, algo = algorithm)

  check.logical(debug)

  # check the cluster.
  cluster = check.cluster(cluster)

  if (!is.null(cluster)) {

    # enter in parallel mode.
    parallel = TRUE
    # set up the slave processes.
    slaves.setup(cluster)
    # disable debugging, the slaves do not cat() here.
    if (debug) {

      warning("disabling debugging output for parallel computing.")
      debug = FALSE

    }#THEN

  }#THEN

  if (algorithm == "direct.lingam") {

    maximize = check.dlingam.maximize(maximize)
    maximize.args = check.dlingam.maximize.args(maximize,
                      extra.args = maximize.args, data = data)

    dag = dlingam.backend(data = data, whitelist = whitelist,
            blacklist = blacklist, mi = mi, maximize = maximize,
            maximize.args = maximize.args, debug = debug, cluster = cluster)

  }#THEN

  dag$learning = list(whitelist = whitelist, blacklist = blacklist, test = mi,
                      algo = algorithm, ordering = dag$learning$ordering,
                      ntests = dag$learning$ntests, args = list())
  # add tests performed by the slaves to the test counter.
  if (parallel)
    dag$learning$ntests = dag$learning$ntests +
      sum(unlist(parallel::clusterEvalQ(cluster, test.counter())))

  invisible(dag)

}#LINGAM.LEARNERS
