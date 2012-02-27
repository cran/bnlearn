
# measure the strength of the arcs in a directed graph.
arc.strength = function(x, data, criterion = NULL, ..., debug = FALSE) {

  # check x's class.
  check.bn(x)
  # arc strength is undefined in partially directed graphs.
  if (is.pdag(x$arcs, names(x$nodes)) && !identical(criterion, "bootstrap"))
    stop("the graph is only partially directed.")
  # check the data are there.
  check.data(data)
  # check the network against the data.
  check.bn.vs.data(x, data)
  # check debug.
  check.logical(debug)
  # check criterion.
  if (is.null(criterion)) {

    # if no criterion is specified use either the default one or the
    # one used by the learning algorithm.
    if (x$learning$test == "none")
      criterion = check.test(criterion, data)
    else
      criterion = x$learning$test

  }#THEN
  else if (identical(criterion, "bootstrap")) {

    # nothing to do, move along.

  }#THEN
  else  {

    criterion = check.criterion(criterion, data)

  }#ELSE

  # set the test/score counter.
  assign(".test.counter", 0, envir = .GlobalEnv)

  # expand and sanitize score-specific arguments and the alpha threshold.
  if (criterion %in% available.tests) {

    # sanitize the alpha threshold.
    alpha = check.alpha(list(...)$alpha, network = x)

    # sanitize B (the number of bootstrap/permutation samples).
    B = check.B(list(...)$B, criterion)

    # warn about unused arguments.
    check.unused.args(list(...), c("alpha", "B"))

    res = arc.strength.test(network = x, data = data, alpha = alpha,
            test = criterion, B = B, debug = debug)

    # add extra information for strength.plot().
    res = structure(res, mode = "test", threshold = alpha)

  }#THEN
  else if (criterion %in% available.scores) {

    # expand and sanitize score-specific arguments.
    extra.args = check.score.args(score = criterion, network = x,
                   data = data, extra.args = list(...))

    res = arc.strength.score(network = x, data = data, score = criterion,
            extra = extra.args, debug = debug)

    # add extra information for strength.plot().
    res = structure(res, mode = "score", threshold = 0)

  }#THEN
  else if (criterion == "bootstrap") {

    # expand and check bootstrap-specific arguments.
    extra.args = check.bootstrap.args(list(...), network = x, data = data)

    res = arc.strength.boot(data = data, R = extra.args$R,
            m = extra.args$m, algorithm = extra.args[["algorithm"]],
            algorithm.args = extra.args[["algorithm.args"]], arcs = NULL,
            cpdag = FALSE, debug = debug)

    # compute the confidence threshold before subsetting.
    t = threshold(res)
    # extract the subset of arcs of interest.
    res = match.arcs.and.strengths(x$arcs, names(x$nodes), res, keep = TRUE)

    # add extra information for strength.plot(), and drop the column
    # with the direction confidence.
    res = structure(res[, 1:3], mode = "bootstrap", threshold = t)

  }#THEN

  return(structure(res, row.names = seq(nrow(res)), class = c("bn.strength", class(res))))

}#ARC.STRENGTH

# compute the strength of all possible arcs from a list of network
# structures/arc sets.
custom.strength = function(networks, nodes, weights = NULL, cpdag = TRUE, debug = FALSE) {

  # check debug.
  check.logical(debug)
  # check cpdag.
  check.logical(cpdag)
  # check the node labels.
  check.nodes(nodes, min.nodes = 3)
  # check networks.
  check.customlist(networks, nodes = nodes)
  # check the weights.
  weights = check.weights(weights, length(networks))

  res = arc.strength.custom(custom.list = networks, nodes, cpdag = cpdag,
          arcs = NULL, weights = weights, debug = debug)

  # add extra information for strength.plot().
  res = structure(res, mode = "bootstrap", threshold = threshold(res),
          class = c("bn.strength", class(res)))

  return(res)

}#CUSTOM.STRENGTH

# build the averaged network structure using arc strengths and a
# significance threshold.
averaged.network = function(strength, nodes, threshold) {

  # check the strength parameter.
  check.bn.strength(strength)
  # this works only with bootstrapped networks.
  if (attributes(strength)$mode != "bootstrap")
    stop("only arc strength computed from bootstrapped networks are supported.")
  # check the strength threshold.
  threshold = check.threshold(threshold, strength)
  # check nodes.
  if (missing(nodes)) {

    # use the bn.strength object to get a node set.
    nodes = unique(c(strength[, "from"], strength[, "to"]))

  }#THEN
  else {

    # sanitize the node set.
    check.nodes(nodes = nodes, min.nodes = 3)
    # double-check whther the bn.strength object agrees with the node set.
    check.bn.strength(strength, nodes = nodes)

  }#ELSE

  avg = averaged.network.backend(strength = strength, nodes = nodes,
          threshold = threshold)

  # add the metadata for the print() method.
  avg$learning$algo = "averaged"
  avg$learning$args = list(threshold = threshold)

  return(avg)

}#AVERAGED.NETWORK

