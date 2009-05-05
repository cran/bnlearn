
# Parameter sanitization for the constraint-based learning algorithms.
bnlearn = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
    test = "mi", alpha = 0.05, B = NULL, method = "gs", debug = FALSE, 
    optimized = TRUE, strict = TRUE, undirected = FALSE) {

  assign(".test.counter", 0, envir = .GlobalEnv)

  res = NULL
  available.methods = c("gs", "iamb", "fast-iamb", "inter-iamb", "mmpc")
  supported.clusters = c("MPIcluster", "PVMcluster","SOCKcluster")
  cluster.aware = FALSE

  # check the data are there.
  check.data(x)
  # check the algorithm.
  if (!(method %in% available.methods))
    stop(paste("valid values for method are:",
           paste(available.methods, collapse = " ")))
  # check test labels.
  test = check.test(test, x)
  # check the logical flags (debug, strict, optimized, undirected, direction).
  check.logical(debug)
  check.logical(strict)
  check.logical(optimized)
  check.logical(undirected)
  # check alpha.
  alpha = check.alpha(alpha)
  # check B (the number of bootstrap/permutation samples).
  B = check.B(B, test)

  # check cluster.
  if (!is.null(cluster)) {

    if (!(any(class(cluster) %in% supported.clusters)))
      stop("cluster is not a valid cluster object.")
    else if (!(require(snow)))
      stop("Can't find required packages: snow")
    else if (!isClusterRunning(cluster))
      stop("the cluster is stopped.")
    else {

      # enter in cluster-aware mode.
      cluster.aware = TRUE
      # set the test counter in all the cluster nodes.
      clusterEvalQ(cluster, assign(".test.counter", 0, envir = .GlobalEnv))
      # disable debugging, the slaves do not cat() here.
      if (debug) {

        warning("disabling debugging output in cluster-aware mode.")
        debug = FALSE

      }#THEN

    }#ELSE

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
  else if (method == "fast-iamb") {

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
  else if (method == "inter-iamb") {

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

  if (undirected) {

    # save the status of the learning algorithm.
    arcs = nbr2arcs(mb)
    learning = list(whitelist = whitelist, blacklist = blacklist, 
      test = test, args = list(alpha = alpha),
      ntests = get(".test.counter", envir = .GlobalEnv))

    # include also the number of permutations/bootstrap samples
    # if it makes sense.
    if (!is.null(B))
      learning$args$B = B

    res = list(learning = learning, 
      nodes = cache.structure(names(mb), arcs = arcs), arcs = arcs)

  }#THEN
  else {

    # recover some of the arc directions.
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

  invisible(structure(res, class = "bn"))

}#BNLEARN

# Parameter sanitization for the score-based learning algorithms.
greedy.search = function(x, start = NULL, whitelist = NULL, blacklist = NULL,
    score = "k2", heuristic = "hc", ..., debug = FALSE, restart = 0,
    perturb = 1, max.iter = Inf, optimized = FALSE) {

  # check the data are there.
  check.data(x)
  # check the score label.
  score = check.score(score, x)
  # check debug.
  check.logical(debug)
  # check restart and perturb.
  check.restart(restart, perturb)
  # check the max.iter parameter
  if ((max.iter != Inf) && !is.positive.integer(max.iter))
    stop("the maximum number of iterations must be a positive integer number.")

  # sanitize whitelist and blacklist, if any.
  whitelist = build.whitelist(whitelist, names(x))
  blacklist = build.blacklist(blacklist, whitelist, names(x))
  # if there is no preseeded network, use an empty one.
  if (is.null(start))
    start = empty.graph(nodes = names(x))
  else {

    # check start's class.
    if (!is(start, "bn"))
      stop("x must be an object of class 'bn'.")
    # set all nodes as updated if the preseed network is not empty,
    # so that all cache lookups are skipped.
    if (nrow(start$arcs) > 0)
      start$updates = array(rep(0, length(start$nodes)),
                        dimnames = list(names(start$nodes)))

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
  if (!is.acyclic.backend(start$arcs, names(start$nodes), directed = TRUE))
    stop("the preseeded graph contains cycles.")

  # expand and sanitize score-specific arguments.
  extra.args = check.score.args(score = score, network = start,
                 data = x, extra.args = list(...))

  # create the test counter in .GlobalEnv.
  assign(".test.counter", 0, envir = .GlobalEnv)

  if (heuristic == "hc") {

    if (optimized) {

      res = hill.climbing.optimized(x = x, start = start,
        whitelist = whitelist, blacklist = blacklist, score = score,
        extra.args = extra.args, restart = restart, perturb = perturb,
        max.iter = max.iter, debug = debug)

    }#THEN
    else {

      res = hill.climbing(x = x, start = start, whitelist = whitelist,
        blacklist = blacklist, score = score, extra.args = extra.args,
        restart = restart, perturb = perturb, max.iter = max.iter,
        debug = debug)

    }#ELSE

  }#THEN

  # set the metadata of the network.
  res$learning$algo = heuristic
  res$learning$ntests = get(".test.counter", envir = .GlobalEnv)
  res$learning$test = score
  res$learning$args = extra.args

  invisible(res)

}#GREEDY.SEARCH

