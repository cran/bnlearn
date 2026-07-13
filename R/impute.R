
# missing data imputation backend.
impute.backend = function(fitted, data, cluster = NULL, method, extra.args,
    restrict.to = NULL, prob = FALSE, debug = FALSE) {

  # nothing to do if there are no observations.
  if (nrow(data) == 0)
    return(data)

  # individual incomplete observations are imputed independently of each other,
  # so they can be processed in parallel.
  if (!is.null(cluster))
    splits = parallel::clusterSplit(cluster, seq(nrow(data)))
  else
    splits = list(seq(nrow(data)))

  if (method == "parents")
    fun = impute.backend.parents
  else if (method == "bayes-lw")
    fun = impute.backend.likelihood.weighting
  else if (method == "exact")
    fun = impute.backend.exact

  imputed = smartSapply(cluster, splits, function(ids) {
    fun(fitted, data = data[ids, , drop = FALSE], extra.args = extra.args,
        restrict.to = restrict.to, prob = prob, debug = debug)
  })

  if (!is.null(cluster)) {

    # rbind() drops attributes; rebuild the prob matrix from the splits first.
    if (prob)
      probabilities = do.call("cbind", lapply(imputed, attr, "prob"))
    imputed = do.call("rbind", imputed)
    if (prob)
      attr(imputed, "prob") = probabilities

  }#THEN
  else
    imputed = imputed[[1]]

  # ensure the metadata are up-to-date (this keeps the "prob" attribute).
  attr(imputed, "metadata") = collect.metadata(imputed)

  return(imputed)

}#IMPUTE.BACKEND

# missing data imputation with maximum likelihood predictions.
impute.backend.parents = function(fitted, data, extra.args, restrict.to = NULL,
    prob = FALSE, debug = FALSE) {

  # check the variables in topological order, to ensure parents are complete.
  for (i in topological.ordering(fitted)) {

    if (debug)
      cat("* checking node", i, ".\n")

    # if there is no missing value, nothing to do.
    missing = is.na(data[, i])

    if (!any(missing))
      next

    if (debug)
      cat("  > found", how.many(missing), "missing value(s).\n")

    # extract the data from the parents of the node.
    predict.from = data[missing, parents(fitted, i), drop = FALSE]

    # call predict.backend() so that arguments are not sanitized again.
    data[missing, i] =
      predict.backend(fitted = fitted, node = i, data = predict.from,
        cluster = NULL, method = "parents", prob = FALSE, debug = FALSE)

  }#FOR

  return(data)

}#IMPUTE.BACKEND.PARENTS

# collapse the per-observation posteriors into a "prob" matrix.
collapse.imputation.probabilities = function(probs, node.levels) {

  k = length(node.levels)

  columns = lapply(probs, function(p) {

    column = structure(rep(NA_real_, k), names = node.levels)
    if (!is.null(p))
      column[names(p)] = p

    return(column)

  })

  matrix(unlist(columns), nrow = k, ncol = length(probs),
    dimnames = list(node.levels, NULL))

}#COLLAPSE.IMPUTATION.PROBABILITIES

# missing data imputation with maximum likelihood predictions.
impute.backend.likelihood.weighting = function(fitted, data, extra.args,
    restrict.from = NULL, restrict.to = NULL, prob = FALSE, debug = FALSE) {

  nodes = names(data)
  n = extra.args$n

  # if there is no missing value, nothing to do.
  missing = !complete.cases(data)

  # posterior probabilities of the imputed values, one vector per observation.
  if (prob)
    probs = vector("list", nrow(data))

  for (j in which(missing)) {

    from = nodes[which(!is.na(data[j, ]))]
    to = setdiff(nodes, from)

    # further reduce the set of nodes to impute from, if required.
    if (!is.null(restrict.from))
      from = intersect(from, restrict.from)

    # further reduce the set of nodes to impute, if required.
    if (!is.null(restrict.to)) {

      to = intersect(to, restrict.to)
      if (length(to) == 0)
        next

    }#THEN

    # use the observed part of the observation as the evidence.
    if (length(from) == 0)
      evidence = TRUE
    else
      evidence = lapply(data[j, from, drop = FALSE],
                   function(x) if (is.factor(x)) as.character(x) else x)

    if (debug) {

      cat("* observation", j, ", imputing", to, "from: \n")
      print(data.frame(evidence), row.names = FALSE)

    }#THEN

    # simulate the particles and the weights using likelihood weighting.
    particles = conditional.distribution(fitted = fitted, nodes = to,
                  evidence = evidence, method = "lw",
                  extra = list(from = from, n = n), debug = FALSE)

    # impute by posterior mode (discrete variables) or posterior expectation
    # (continuous variables); discard missing weights.
    w = attr(particles, "weights")
    w[is.na(w)] = 0

    estimates = lapply(particles, function(x, w) {

      if (is.factor(x))
        names(which.max(by(w, INDICES = x, FUN = sum)))
      else if (is.numeric(x))
        weighted.mean(x = x, w = w)

    }, w = w)

    # make sure unsuccessful imputations can be handled in the assignment later.
    failures = sapply(estimates, function(x) is.na(x) || is.null(x))
    if (any(failures))
      estimates[failures] = NA

    # posterior of the single discrete target, from its weighted level counts.
    if (prob && (length(to) == 1) && is.factor(particles[[to]]) &&
        !failures[[to]]) {

      counts = tapply(w, particles[[to]], sum)[levels(particles[[to]])]
      counts[is.na(counts)] = 0
      probs[[j]] = counts / sum(counts)

    }#THEN

    if (debug) {

      cat("  > imputed value:", "\n")
      print(data.frame(estimates), row.names = FALSE)

    }#THEN

    data[j, to] = estimates

  }#FOR

  # collapse the per-observation probabilities into the "prob" matrix.
  if (prob) {

    target.levels = dimnames(fitted[[restrict.to]]$prob)[[1]]
    attr(data, "prob") =
      collapse.imputation.probabilities(probs, target.levels)

  }#THEN

  return(data)

}#IMPUTE.BACKEND.LIKELIHOOD.WEIGHTING

# missing data imputation with exact inference.
impute.backend.exact = function(fitted, data, extra.args, restrict.from = NULL,
    restrict.to = NULL, prob = FALSE, debug = FALSE) {

  if (is(fitted, c("bn.fit.dnet", "bn.fit.onet", "bn.fit.donet"))) {

    exact.discrete.imputation(fitted = fitted, data = data,
      restrict.from = restrict.from, restrict.to = restrict.to, prob = prob,
      debug = debug)

  }#THEN
  else if (is(fitted, "bn.fit.gnet")) {

    exact.gaussian.imputation(fitted = fitted, data = data,
      restrict.from = restrict.from, restrict.to = restrict.to, prob = prob,
      debug = debug)

  }#THEN
  else if (is(fitted, "bn.fit.cgnet")) {

    stop("'bn.fit.cgnet' networks are not supported.")

  }#ELSE

}#IMPUTE.BACKEND.EXACT

# missing data imputation with junction trees and belief propagation.
exact.discrete.imputation = function(fitted, data, restrict.from = NULL,
    restrict.to = NULL, prob = FALSE, debug = FALSE) {

  nodes = names(data)

  # check whether gRain is loaded.
  check.and.load.package("gRain")

  jtree = from.bn.fit.to.grain(fitted, compile = TRUE)

  # for each incomplete observation...
  missing = !complete.cases(data)

  # posterior probabilities of the imputed values, one vector per observation.
  if (prob)
    probs = vector("list", nrow(data))

  for (j in which(missing)) {

    # ... separate observed and unobserved values...
    from = nodes[which(!is.na(data[j, ]))]
    to = setdiff(nodes, from)

    # further reduce the set of nodes to impute from, if required.
    if (!is.null(restrict.from))
      from = intersect(from, restrict.from)

    # further reduce the set of nodes to impute, if required.
    if (!is.null(restrict.to)) {

      to = intersect(to, restrict.to)
      if (length(to) == 0)
        next

    }#THEN

    if (debug) {

      cat("* observation", j, ", imputing", to, "from: \n")
      print(data[j, from, drop = FALSE], row.names = FALSE)

    }#THEN

    # ... split the imputation into more manageable chunks if possible...
    queries = query.partitioning(fitted, event = to, evidence = from,
                debug = debug)

    for (q in queries) {

      # ... impute the missing values with their most probable explanation...
      imputed = mpe.discrete.query(jtree, event = q$event,
                  evidence = q$evidence, value = data[j, , drop = FALSE])

      # ... and save the imputed observation back into the data set, if it was
      # possible to compute them at all.
      if (!is.null(imputed)) {

        data[j, q$event] = imputed

        # posterior of the single target, reusing the MPE's joint table.
        if (prob && !is.null(restrict.to) && (length(q$event) == 1) &&
            (q$event == restrict.to)) {

          posterior = attr(imputed, "prob")
          probs[[j]] = structure(as.numeric(posterior),
                         names = dimnames(posterior)[[restrict.to]])

        }#THEN

      }#THEN

    }#FOR

    if (debug) {

      cat("  > imputed value:", "\n")
      print(data.frame(data[j, to, drop = FALSE]), row.names = FALSE)

    }#THEN

  }#FOR

  # collapse the per-observation probabilities into the "prob" matrix.
  if (prob) {

    target.levels = dimnames(fitted[[restrict.to]]$prob)[[1]]
    attr(data, "prob") =
      collapse.imputation.probabilities(probs, target.levels)

  }#THEN

  return(data)

}#EXACT.DISCRETE.IMPUTATION

# missing data imputation with closed-form Gaussian results.
exact.gaussian.imputation = function(fitted, data, restrict.from = NULL,
    restrict.to = NULL, prob = FALSE, debug = FALSE) {

  nodes = names(data)

  # get the global distribution...
  mvn = gbn2mvnorm.backend(fitted)
  # ... and check that it is not singular.
  if (anyNA(mvn$mu) || anyNA(mvn$sigma))
    return(data)

  missing = !complete.cases(data)

  # for each incomplete observation...
  for (j in which(missing)) {

    # ... separate observed and unobserved values...
    from = nodes[which(!is.na(data[j, ]))]
    to = setdiff(nodes, from)

    # further reduce the set of nodes to impute from, if required.
    if (!is.null(restrict.from))
      from = intersect(from, restrict.from)

    # further reduce the set of nodes to impute, if required.
    if (!is.null(restrict.to)) {

      to = intersect(to, restrict.to)
      if (length(to) == 0)
        next

    }#THEN

    if (debug) {

      cat("* observation", j, ", imputing", to, "from: \n")
      print(data[j, from, drop = FALSE], row.names = FALSE)

    }#THEN

    # ... split the imputation into more manageable chunks if possible...
    queries = query.partitioning(fitted, event = to, evidence = from,
                debug = debug)

    for (q in queries) {

      # ...and impute the missing values with their most probable explanation.
      data[j, q$event] =
        mpe.gaussian.query(mvn, event = q$event, evidence = q$evidence,
                        value = data[j, q$evidence, drop = FALSE])

    }#THEN

    if (debug) {

      cat("  > imputed value:", "\n")
      print(data[j, to, drop = FALSE], row.names = FALSE)

    }#THEN

  }#FOR

  return(data)

}#EXACT.GAUSSIAN.IMPUTATION
