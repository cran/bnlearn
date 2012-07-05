
# simple parametric and nonparametric bootstrap implementation.
bootstrap.backend = function(data, statistic, R, m, sim = "ordinary",
    algorithm, algorithm.args = list(), statistic.args = list(),
    cluster = NULL, debug = FALSE) {

  # allocate the result list.
  res = as.list(seq(R))
  # allocate the bayesian network to use for parametric bootstrap.
  net = NULL

  # initialize the bayesian network used by the paramentric bootstrap.
  if (sim == "parametric") {

    # the mbde score requires knowledge of which observations have been
    # manipulated, and this is clearly not possible when samples are
    # generated.
    if (algorithm.args$score == "mbde")
      stop("the 'mbde' score is incompatible with parametric bootstrap.")

    if (algorithm %in% always.dag.result) {

      net = do.call(algorithm, c(list(x = data), algorithm.args))

      if (debug) {

        cat("* initial network for parametric bootstrap is:\n")
        print(net)

      }#THEN

    }#THEN
    else {

      stop(paste("this learning algorithm may result in a partially",
             "directed dag, which is not handled by parametric bootstrap."))

    }#ELSE

  }#THEN

  bootstrap.step = function(r, data, m, net, sim, algorithm, algorithm.args,
      statistic, statistic.args, debug) {

    if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("* bootstrap replicate", r, ".\n")

    }#THEN

    # generate the r-th bootstrap sample.
    if (sim == "ordinary") {

      # perform the resampling with replacement.
      resampling = sample(nrow(data), m, replace = TRUE)

      # user-provided lists of manipulated observations for the mbde score must
      # be remapped to match the bootstrap sample.
      if ((algorithm.args$score == "mbde") && !is.null(algorithm.args$exp)) {

        algorithm.args$exp = lapply(algorithm.args$exp, function(x) {

          x = match(x, resampling)
          x = x[!is.na(x)]

        })

      }#THEN

      # generate the bootstrap sample.
      replicate = data[resampling, , drop = FALSE]

    }#THEN
    else if (sim == "parametric") {

      if (is.data.discrete(data))
        replicate = rbn.discrete(x = net, n = m, data = data)
      else
        replicate = rbn.continuous(x = net, n = m, data = data)

    }#ELSE

    if (debug)
      cat("* learning bayesian network structure.\n")

    # learn the network structure from the bootstrap sample.
    bn = do.call(algorithm, c(list(x = replicate), algorithm.args))

    if (debug) {

      print(bn)
      cat("* applying user-defined statistic.\n")

    }#THEN

    # apply the user-defined function to the newly-learned bayesian network;
    # the bayesian network is passed as the first argument hoping it will end
    # at the right place thanks to the positional matching.
    res = do.call(statistic, c(list(bn), statistic.args))

    if (debug) {

      cat("  > the function returned:\n")
      print(res)

    }#THEN

    return(res)

  }#BOOTSTRAP.STEP

  if (!is.null(cluster)) {

    res = parLapply(cluster, res, bootstrap.step, data = data, m = m, net = net,
      sim = sim, algorithm = algorithm, algorithm.args = algorithm.args,
      statistic = statistic, statistic.args = statistic.args, debug = debug)

  }#THEN
  else {

    res = lapply(res, bootstrap.step, data = data, m = m, net = net, sim = sim,
      algorithm = algorithm, algorithm.args = algorithm.args, statistic = statistic,
      statistic.args = statistic.args, debug = debug)

  }#ELSE

  return(res)

}#BOOTSTRAP.BACKEND

# model averaging for bootstrapped network structures.
averaged.network.backend = function(strength, nodes, threshold) {

  e = empty.graph(nodes)

  # arcs with a strength of one should always be selected, regardless of
  # the threshold.
  significant = (strength$strength > threshold) | (strength$strength == 1)

  # filter also the direction if present in the bn.strength object.
  if ("direction" %in% names(strength))
    significant = significant & (strength$direction >= 0.5)

  # nothing to see, move along.
  if (!any(significant))
    return(e)

  candidate.arcs = as.matrix(strength[significant, c("from", "to"), drop = FALSE])

  if (all(which.undirected(candidate.arcs))) {

    # update the arcs of the network, no cycles.
    e$arcs = candidate.arcs

  }#THEN
  else {

    # update the arcs of the network, minding cycles.
    e$arcs = .Call("smart_network_averaging",
                   arcs = candidate.arcs,
                   nodes = nodes,
                   weights = strength$strength[significant],
                   PACKAGE = "bnlearn")

  }#ELSE

  # update the network structure.
  e$nodes = cache.structure(nodes, arcs = candidate.arcs)

  return(e)

}#AVERAGED.NETWORK.BACKEND

