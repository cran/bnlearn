
# backend for sampling from conditional probability distributions.
conditional.distribution = function(fitted, nodes, evidence, method, extra,
    cluster = NULL, debug = FALSE) {

  # consider only the upper closure of event and evidence to reduce the number
  # of variables in the Monte Carlo simulation.
  fitted = reduce.fitted(fitted = fitted, event = nodes, evidence = evidence,
             nodes = extra$query.nodes, method = method, debug = debug)

  if (method == "ls")
    distribution = logic.distribution
  else if (method == "lw")
    distribution = weighting.distribution

  if (!is.null(cluster)) {

    # get the number of slaves.
    s = nSlaves(cluster)
    # divide the number of particles among the slaves.
    batch = n = ceiling(extra$n / s)

    results = parallel::parLapplyLB(cluster, seq(s),
      function(x) {

        distribution(fitted = fitted, nodes = nodes, evidence = evidence,
          n = n, batch = batch, debug = debug)

      })

    return(do.call(rbind, results))

  }#THEN
  else {

    distribution(fitted = fitted, nodes = nodes, evidence = evidence,
      n = extra$n, batch = extra$batch, debug = debug)

  }#ELSE

}#CONDITIONAL.DISTRIBUTION

# backend for conditional probability queries.
conditional.probability.query = function(fitted, event, evidence, method,
    extra, cluster = NULL, debug = FALSE) {

  # consider only the upper closure of event and evidence to reduce the number
  # of variables in the Monte Carlo simulation.
  fitted = reduce.fitted(fitted = fitted, event = event, evidence = evidence,
             nodes = extra$query.nodes, method = method, debug = debug)

  if (method == "ls")
    sampling = logic.sampling
  else if (method == "lw")
    sampling = weighting.sampling

  if (!is.null(cluster)) {

    # get the number of slaves.
    s = nSlaves(cluster)
    # divide the number of particles among the slaves.
    batch = n = ceiling(extra$n / s)

    results = parallel::parSapplyLB(cluster, seq(s),
      function(x) {

        sampling(fitted = fitted, event = event, evidence = evidence,
         n = n, batch = batch, debug = debug)

      })

    return(mean(results))

  }#THEN
  else {

    sampling(fitted = fitted, event = event, evidence = evidence,
      n = extra$n, batch = extra$batch, debug = debug)

  }#ELSE

}#CONDITIONAL.PROBABILITY.QUERY

# create an empty data frame from a bn.fit object.
fit.dummy.df = function(fitted, nodes) {

  dummy = sapply(nodes, function(x) {

    node = fitted[[x]]

    if (is(node, "bn.fit.dnode"))
      return(factor(character(0), levels = dimnames(node$prob)[[1]]))
    else if (is(node, "bn.fit.onode"))
      return(ordered(character(0), levels = dimnames(node$prob)[[1]]))
    else if (is(node, "bn.fit.gnode"))
      return(numeric(0))

  })

  return(.data.frame(dummy))

}#FIT.DUMMY.DF

# reduce a bn.fit object to the upper closure of event and evidence nodes.
reduce.fitted = function(fitted, event, evidence, nodes, method, debug = FALSE) {

  if (is.null(nodes)) {

    # find out which nodes are involved in the event and the evidence and
    # construct their upper closure.
    nodes = names(fitted)
    nodes.event = nodes[nodes %in% explode(event)]
    nodes.evidence = nodes[nodes %in% explode(evidence)]
    upper.closure =
      topological.ordering(fitted, start = union(nodes.evidence, nodes.event),
                           reverse = TRUE)

    # check whether something went horribly wrong in getting the node labels
    # and the upper closure.
    if (((length(upper.closure) == 0) || any(upper.closure %!in% nodes)) ||
        (!identical(evidence, TRUE) && (length(nodes.evidence) == 0)) ||
        (!identical(event, TRUE) && (length(nodes.event) == 0))) {

      # use all the nodes when in doubt.
      event = evidence = TRUE
      upper.closure = nodes

      if (debug) {

        cat("* checking which nodes are needed.\n")
        cat("  > unable to identify query nodes.\n")

      }#THEN

    }#THEN
    else {

      if (debug) {

        cat("* checking which nodes are needed.\n")
        cat("  > event involves the following nodes:", nodes.event, "\n")
        cat("  > evidence involves the following nodes:", nodes.evidence, "\n")
        cat("  > upper closure is '", upper.closure, "'\n")

      }#THEN

    }#ELSE

  }#THEN
  else {

    # construct the upper closure of the query nodes.
    upper.closure = topological.ordering(fitted, start = nodes, reverse = TRUE)

    if (debug) {

      cat("* using specified query nodes.\n")
      cat("  > upper closure is '", upper.closure, "'\n")

    }#THEN

  }#ELSE

  # check whether the upper closure is correct: tricky expressions are not
  # always handled correctly by explode().
  dummy = fit.dummy.df(fitted, upper.closure)
  # testing evidence (TRUE or an expression in LS, a list in LW).
  if (!(is.language(evidence) || identical(evidence, TRUE)))
    try.evidence = TRUE
  else
    try.evidence = try(eval(evidence, dummy), silent = TRUE)
  # testing event (it's the label nodes in cpdist; TRUE or an expression
  # in cpquery).
  if (!(is.language(event) || identical(event, TRUE)))
    try.event = TRUE
  else
    try.event = try(eval(event, dummy), silent = TRUE)

  # create the subgraph corresponding to the upper closure, keeping the ordering
  # of the node in the bn.fit object the same to avoid potentially affecting
  # random simulations later.
  if (is.logical(try.event) && is.logical(try.evidence)) {

    if (debug)
      cat("  > generating observations from", length(upper.closure), "/",
        length(fitted), "nodes.\n")

    fitted =
      structure(fitted[names(fitted) %in% upper.closure], class = class(fitted))

  }#THEN
  else {

    if (debug)
      cat("  > unable use the upper closure, using the whole network.\n")

  }#ELSE

  return(fitted)

}#REDUCE.FITTED

# compute conditional probabilities with forward/logic sampling.
logic.sampling = function(fitted, event, evidence, n, batch, debug = FALSE) {

  cpxe = cpe = 0L
  filtered = logical(n)
  matching = logical(n)
  r = logical(n)

  # count how many complete batches we have to generate.
  nbatches = n %/% batch
  # count how many observations are in the last one.
  last.one = n %% batch

  for (m in c(rep(batch, nbatches), last.one)) {

    # do a hard reset of generated.data, so that the memory used by the data
    # set generated in the previous iteration can be garbage-collected and
    # reused _before_ rbn() returns.
    generated.data = NULL

    # generate random data from the bayesian network.
    if (m > 0) {

      generated.data = rbn.backend(x = fitted, n = m)

    }#THEN
    else
      break

    if (debug)
      cat("* generated", m, "samples from the bayesian network.\n")

    # evaluate the expression defining the evidence.
    if (identical(evidence, TRUE))
      r = rep(TRUE, m)
    else
      r = eval(evidence, generated.data, parent.frame())
    # double-check that this is a logical vector.
    if (!is.logical(r))
      stop("evidence must evaluate to a logical vector.")
    # double-check that it has the right length.
    if (length(r) != m)
      stop("logical vector for evidence is of length ", length(r),
        " instead of ", m, ".")
    # filter out the samples not matching the evidence we assume, and those for
    # which the evidence evaluates to NA because the network is singular.
    filtered = r & !is.na(r)

    # evaluate the expression defining the event.
    if (identical(event, TRUE))
      r = rep(TRUE, m)
    else
      r = eval(event, generated.data, parent.frame())
    # double-check that this is a logical vector.
    if (!is.logical(r))
      stop("event must evaluate to a logical vector.")
    # double-check that it has the right length.
    if (length(r) != m)
      stop("logical vector for event is of length ", length(r),
        " instead of ", m, ".")
    # filter out the samples not matching the event we are looking for, and
    # those for which the event evaluates to NA because the network is singular.
    matching = filtered & r & !is.na(r)
    filtered = filtered & !is.na(r)

    # update the global counters.
    cpe = cpe + how.many(filtered)
    cpxe = cpxe + how.many(matching)

    if (debug) {

      filter.count = how.many(filtered)
      if (!identical(evidence, TRUE))
        cat("  > evidence matches ", filter.count, " samples out of ", m,
          " (p = ", filter.count / m, ").\n", sep = "")
      else
        cat("  > evidence matches ", m, " samples out of ", m,
          " (p = 1).\n", sep = "")

      match.count = how.many(matching)
      fmratio = ifelse(filter.count == 0, 0, match.count / filter.count)
      if (!identical(event, TRUE))
        cat("  > event matches ", match.count, " samples out of ",
          filter.count, " (p = ", fmratio, ").\n", sep = "")
      else
        cat("  > event matches ", filter.count, " samples out of ",
          filter.count, " (p = 1).\n", sep = "")

    }#THEN

  }#FOR

  # prevent divide-by-zero errors.
  result = ifelse(cpe == 0, 0, cpxe / cpe)

  if (debug && (nbatches > 1)) {

    cat("* generated a grand total of", n, "samples.\n")
    cat("  > event matches ", cpxe, " samples out of ", cpe,
      " (p = ", result, ").\n", sep = "")

  }#THEN

  return(result)

}#LOGIC.SAMPLING

# generate random observations from conditional distributions with forward/logic
# sampling.
logic.distribution = function(fitted, nodes, evidence, n, batch, debug = FALSE) {

  filtered = logical(n)
  result = NULL

  # count how many complete batches we have to generate.
  nbatches = n %/% batch
  # count how many observations are in the last one.
  last.one = n %% batch

  for (m in c(rep(batch, nbatches), last.one)) {

    # do a hard reset of generated.data, so that the memory used by the data
    # set generated in the previous iteration can be garbage-collected and
    # reused _before_ rbn() returns.
    generated.data = NULL

    # generate random data from the bayesian network.
    if (m > 0)
      generated.data = rbn.backend(x = fitted, n = m)
    else
      break

    if (debug)
      cat("* generated", m, "samples from the bayesian network.\n")

    # evaluate the expression defining the evidence.
    r = eval(evidence, generated.data, parent.frame())
    # double-check that this is a logical vector.
    if (!is.logical(r))
      stop("evidence must evaluate to a logical vector.")
    # double-check that it has the right length.
    if ((length(r) != 1) && (length(r) != m))
      stop("logical vector for evidence is of length ", length(r),
        " instead of ", m, ".")
    # filter out the samples not matching the evidence we assume.
    filtered = r & !is.na(r)

    if (debug) {

      filter.count = how.many(filtered)
      if (!identical(evidence, TRUE))
        cat("  > evidence matches ", filter.count, " samples out of ", m,
          " (p = ", filter.count / m, ").\n", sep = "")
      else
        cat("  > evidence matches ", m, " samples out of ", m,
          " (p = 1).\n", sep = "")

    }#THEN

    # update the return value.
    result = rbind(result, generated.data[filtered, nodes, drop = FALSE])

  }#FOR

  # reset the row names.
  rownames(result) = NULL
  # set attributes for later use.
  class(result) = c("bn.cpdist", class(result))
  attr(result, "method") = "ls"

  if (debug && (nbatches > 1))
    cat("* generated a grand total of", n, "samples.\n")

  return(result)

}#LOGIC.DISTRIBUTION

# compute conditional probabilities with likelihood weighting.
weighting.sampling = function(fitted, event, evidence, n, batch, debug = FALSE) {

  cpxe = cpe = 0
  matching = logical(n)
  r = logical(n)

  # count how many complete batches we have to generate.
  nbatches = n %/% batch
  # count how many observations are in the last one.
  last.one = n %% batch

  weights = function(data) {

    if (isTRUE(evidence))
      return(rep(1, nrow(data)))
    else
      .Call(call_lw_weights,
            fitted = fitted,
            data = data,
            keep = names(evidence),
            debug = debug)

  }#WEIGHTS

  for (m in c(rep(batch, nbatches), last.one)) {

    # do a hard reset of generated.data, so that the memory used by the data
    # set generated in the previous iteration can be garbage-collected and
    # reused _before_ rbn() returns.
    generated.data = NULL

    # generate random data from the bayesian network.
    if (m > 0)
      generated.data = rbn.backend(x = fitted, fix = evidence, n = m)
    else
      break

    if (debug)
      cat("* generated", m, "samples from the bayesian network.\n")

    # evaluate the expression defining the event.
    if (identical(event, TRUE))
      r = rep(TRUE, m)
    else
      r = eval(event, generated.data, parent.frame())
    # double-check that this is a logical vector.
    if (!is.logical(r))
      stop("event must evaluate to a logical vector.")
    # double-check that it has the right length.
    if (length(r) != m)
      stop("logical vector for event is of length ", length(r),
        " instead of ", m, ".")
    # filter out the samples not matching the event we are looking for.
    matching = r & !is.na(r)

    # compute the probabilities and use them as weigths.
    attr(generated.data, "metadata") = collect.metadata(generated.data)
    w = weights(generated.data)
    cpe = cpe + sum(w[!is.na(r)])
    cpxe = cpxe + sum(w[matching])

    if (debug)
      cat("  > event has a probability mass of ", cpxe, " out of ", cpe, ".\n", sep = "")

  }#FOR

  # compute the conditional probability.
  result = cpxe / cpe

  if (debug && (nbatches > 1)) {

    cat("* generated a grand total of", n, "samples.\n")
    cat("  > event has a probability mass of ", cpxe, " out of ", cpe,
        " (p = ", result, ").\n", sep = "")

  }#THEN

  return(result)

}#WEIGHTING.SAMPLING

# generate random observations from conditional distributions with likelihood
# weighting.
weighting.distribution = function(fitted, nodes, evidence, n, batch,
    debug = FALSE) {

  .Call(call_cpdist_lw,
        fitted = fitted,
        nodes = nodes,
        n = as.integer(n),
        fix = evidence,
        debug = debug)

}#WEIGHTING.DISTRIBUTION
