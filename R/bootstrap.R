
# simple nonparametric bootstrap implementation.
bootstrap.backend = function(data, statistic, R, m, algorithm,
    algorithm.args = list(), statistic.args = list(), cluster = NULL,
    debug = FALSE) {

  # allocate the result list.
  res = as.list(seq(R))
  # allocate the bayesian network to use for parametric bootstrap.
  net = NULL
  # check the data early on.
  data.info = check.data(data)

  bootstrap.replicate = function(r, data, m, net, algorithm, algorithm.args,
      statistic, statistic.args, debug) {

    if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("* bootstrap replicate", r, ".\n")

    }#THEN

    # generate the r-th bootstrap sample by resampling with replacement.
    resampling = sample(nrow(data), m, replace = TRUE)

    # user-provided lists of manipulated observations for the mbde score must
    # be remapped to match the bootstrap sample.
    if (!is.null(algorithm.args$score) && (algorithm.args$score == "mbde") &&
          !is.null(algorithm.args$exp)) {

      algorithm.args$exp = lapply(algorithm.args$exp, function(x) {

        x = match(x, resampling)
        x = x[!is.na(x)]

      })

    }#THEN

    # generate the bootstrap sample.
    replicate = data[resampling, , drop = FALSE]

    if (debug)
      cat("* learning bayesian network structure.\n")

    # learn the network structure from the bootstrap sample.
    bn = do.call(algorithm, c(list(x = replicate), algorithm.args))

    if (debug) {

      print(bn)
      cat("* computing user-defined statistic.\n")

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

  }#BOOTSTRAP.REPLICATE

  if (!is.null(cluster)) {

    res = parallel::parLapplyLB(cluster, res, bootstrap.replicate, data = data,
            m = m, net = net, algorithm = algorithm,
            algorithm.args = algorithm.args, statistic = statistic,
            statistic.args = statistic.args, debug = debug)

  }#THEN
  else {

    res = lapply(res, bootstrap.replicate, data = data, m = m, net = net,
            algorithm = algorithm, algorithm.args = algorithm.args,
            statistic = statistic, statistic.args = statistic.args,
            debug = debug)

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
    e$arcs = .Call(call_smart_network_averaging,
                   arcs = candidate.arcs,
                   nodes = nodes,
                   weights = strength$strength[significant])

  }#ELSE

  # update the network structure.
  e$nodes = cache.structure(nodes, arcs = e$arcs)
  # add back illegal arcs, so that cpdag() works correctly.
  if ("illegal" %in% names(attributes(strength)))
    e$learning$illegal = attr(strength, "illegal")

  return(e)

}#AVERAGED.NETWORK.BACKEND

