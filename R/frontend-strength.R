
# measure the strength of the arcs in a directed graph.
arc.strength = function(x, data, criterion = NULL, ..., debug = FALSE) {

  # check x's class and get the node labels.
  check.bn(x)
  nodes = names(x$nodes)
  # arc strength is undefined in partially directed graphs.
  if (is.pdag(x$arcs, nodes))
    stop("the graph is only partially directed.")
  # check the data are there.
  data = check.data(data, allow.missing = TRUE)
  # check the network against the data.
  check.bn.vs.data(x, data)
  # check debug.
  check.logical(debug)
  # check the score or the test or whatever is used to measure strength.
  if (is.null(criterion)) {

    # if no criterion is specified use either the default one or the
    # one used by the learning algorithm.
    if (x$learning$test == "none")
      criterion = check.test(criterion, data = data)
    else
      criterion = x$learning$test

  }#THEN
  else {

    criterion = check.criterion(criterion, data = data)

  }#ELSE

  # check whether the network is valid for the method.
  check.arcs.against.assumptions(x$arcs, data, criterion)

  # set the test/score counter.
  reset.test.counter()

  # expand and sanitize criterion-specific arguments.
  extra.args = list(...)

  if (criterion %in% available.tests) {

    # sanitize the alpha threshold.
    alpha = check.alpha(extra.args$alpha, network = x)
    # check the optional arguments to the test.
    extra.args = check.test.args(test = criterion, extra.args = extra.args,
                   network = x, data = data)

    res = arc.strength.test(network = x, data = data, alpha = alpha,
            test = criterion, extra.args = extra.args, debug = debug)

    # add extra information for strength.plot().
    res = structure(res, nodes = nodes, method = "test", threshold = alpha)

  }#THEN
  else if (criterion %in% available.scores) {

    # expand and sanitize score-specific arguments.
    extra.args = check.score.args(score = criterion, network = x,
                   data = data, extra.args = extra.args, learning = FALSE)

    res = arc.strength.score(network = x, data = data, score = criterion,
            extra = extra.args, debug = debug)

    # add extra information for strength.plot().
    res = structure(res, nodes = nodes, method = "score", threshold = 0)

  }#THEN

  # set the class of the return value.
  res = structure(res, class = c("bn.strength", class(res)))

  # reset the row names if there are rows.
  if (nrow(res) > 0)
    res  = structure(res, row.names = seq(nrow(res)))

  return(res)

}#ARC.STRENGTH

# compute an approximation of arc and direction strength from the Bayes factors
# that can be computed from a single MAP network.
bf.strength = function(x, data, score, ..., debug = FALSE) {

  # check whether Rmpfr is loaded.
  check.and.load.package("Rmpfr")
  # check x's class and get the node labels.
  check.bn(x)
  nodes = names(x$nodes)
  # arc strength is undefined in partially directed graphs.
  if (is.pdag(x$arcs, nodes))
    stop("the graph is only partially directed.")
  # check the data are there.
  data = check.data(data, allow.missing = TRUE)
  # check the network against the data.
  check.bn.vs.data(x, data)
  # check debug.
  check.logical(debug)

  # make sure the score function is suitable for computing a Bayes factor.
  data.type = attr(data, "metadata")$type

  if (missing(score)) {

    if (data.type %in% discrete.data.types)
      score = "bde"
    else if (data.type %in% continuous.data.types)
      score = "bge"
    else if (data.type %in% mixed.data.types)
      score = "bic-cg"

  }#THEN
  else {

    score = check.score(score, data = data,
              allowed = c(available.discrete.bayesian.scores,
                          available.continuous.bayesian.scores,
                          available.omnibus.scores,
                          grep("bic", available.scores, value = TRUE)))

  }#ELSE

  # check whether the network is valid for the score.
  check.arcs.against.assumptions(x$arcs, data, score)

  # expand and sanitize score-specific arguments.
  extra.args = check.score.args(score = score, network = x,
                 data = data, extra.args = list(...), learning = FALSE)

  res = bf.strength.backend(x = x, data = data, score = score, debug = debug,
          extra.args = extra.args)

  # add extra information for strength.plot() and averaged.network().
  res = structure(res, nodes = nodes, method = "bayes-factor",
          threshold = threshold(res), class = c("bn.strength", class(res)))
  attr(res, "illegal") = list.illegal.arcs(nodes = names(data), data = data,
                           criterion = score)

  return(res)

}#BF.STRENGTH

# compute the strength of all possible arcs from a list of network
# structures/arc sets.
custom.strength = function(networks, nodes, weights = NULL, cpdag = TRUE,
    debug = FALSE) {

  illegal = NULL

  # check debug.
  check.logical(debug)
  # check cpdag.
  check.logical(cpdag)
  # check networks.
  if (is(networks, c("bn.kcv", "bn.kcv.list"))) {

    extract = function(x) {

      # regenerate the network structure...
      net = bn.net(x$fitted)
      # ... and restore the information from structure learning, including
      # whitelists and blacklists that are needed later in cpdag.backend().
      net$learning = x$learning

      return(net)

    }#EXTRACT

    if (is(networks, "bn.kcv"))
      networks = lapply(networks, extract)
    else if (is(networks, "bn.kcv.list")) {

      networks = lapply(networks, function(x) lapply(x, extract))
      networks = do.call("c", networks)

    }#THEN

    # extract the node labels from the networks, and ignore the argument.
    if (!missing(nodes))
      warning("the labels from the 'nodes' argument will be ignored.")

    nodes = names(networks[[1]]$nodes)

    # extract the illegal arcs, which are the same for all networks.
    illegal = networks[[1]]$learning$illegal

  }#THEN
  else {

    # check the node labels.
    check.nodes(nodes)
    # check the networks.
    check.customlist(networks, nodes = nodes)

    # extract the illegal arcs, but disregard them if they are not the same for
    # all networks.
    all.illegal = lapply(networks, function(x) {

      if(is(x, "bn"))
        x$learning$illegal
      else
        NULL

    })

    if (all(sapply(all.illegal, all.equal, all.illegal[[1]])))
      illegal = all.illegal[[1]]

  }#ELSE
  # check the weights.
  weights = check.weights(weights, length(networks))

  res = arc.strength.custom(custom.list = networks, nodes = nodes, cpdag = cpdag,
          arcs = NULL, weights = weights, illegal = illegal, debug = debug)

  # add extra information for strength.plot() and averaged.network().
  res = structure(res, nodes = nodes, method = "bootstrap",
          threshold = threshold(res), class = c("bn.strength", class(res)))
  if (!is.null(illegal))
    attr(res, "illegal") = illegal

  return(res)

}#CUSTOM.STRENGTH

# average multiple bn.strength objects.
mean.bn.strength = function(x, ..., weights = NULL) {

  # check the bn.strength objects.
  strength = c(list(x), list(...))

  method = character(length(strength))
  nodes = unique(c(x[, "from"], x[, "to"]))

  for (s in seq_along(strength)) {

    # check the nodes are the same, and that the object has the right structure.
    check.bn.strength(strength[[s]])
    check.bn.strength.vs.bn(strength[[s]], nodes)
    # check the objects have the same number of arcs.
    if (length(strength[[s]][, "from"]) != length(x[, "from"]))
      stop("the bn.strength objects have different numbers of arcs.")
    # check that the objects have the same arcs.
    if (!all(which.listed(as.matrix(strength[[s]][, c("from", "to")]),
                         as.matrix(x[, c("from", "to")]))))
      stop("the bn.strength objects have different arcs.")

    method[s] = attr(strength[[s]], "method")

  }#FOR

  # check all strengths are either score differences, p-values or
  # strength/direction probabilities.
  if (any(method != method[1]))
    stop("all the bn.strength objects must be computed from either a score, ",
      "a conditional independence test, or bootstrap resampling.")

  # check the weights.
  weights = check.weights(weights, length(strength))
  # average the objects.
  res = average.strength(strength, nodes, weights)
  # set the attributes of the return value.
  attributes(res) = attributes(x)
  if (attr(res, "method") == "bootstrap")
    attr(res, "threshold") = threshold(res)

  return(res)

}#MEAN.BN.STRENGTH
