
# missing data imputation backend.
impute.backend = function(fitted, data, method, extra.args, debug) {

  if (method == "parents") {

    impute.backend.parents(fitted = fitted, data = data, debug = debug)

  }#THEN
  else if (method == "bayes-lw") {

    impute.backend.map(fitted = fitted, data = data, n = extra.args$n,
      debug = debug)

  }#THEN

}#IMPUTE.BACKEND

# missing data imputation with maximum likelihood predictions.
impute.backend.parents = function(fitted, data, debug = FALSE) {

  # check the variables in topological order, to ensure parents are complete.
  for (i in topological.ordering(fitted)) {

    if (debug)
      cat("* checking node", i, ".\n")

    # if there is no missing value, nothing to do.
    missing = is.na(data[, i])

    if (!any(missing))
      next

    if (debug)
      cat("  > found", length(which(missing)), "missing value(s).\n")

    # extract the data from the parents of the node.
    predict.from = data[missing, parents(fitted, i), drop = FALSE]

    # call predict() backends so that arguments are not sanitized again.
    if (is(fitted, c("bn.fit.dnet", "bn.fit.onet", "bn.fit.donet"))) {

      data[missing, i] =
        discrete.prediction(node = i, fitted = fitted, data = predict.from,
          prob = FALSE, debug = FALSE)

    }#THEN
    else if (is(fitted, "bn.fit.gnet")) {

      data[missing, i] =
        gaussian.prediction(node = i, fitted = fitted, data = predict.from,
          debug = FALSE)

    }#THEN
    else if (is(fitted, "bn.fit.cgnet")) {

      data[missing, i] =
        mixedcg.prediction(node = i, fitted = fitted, data = predict.from,
          debug = FALSE)

    }#ELSE

  }#FOR

  return(data)

}#IMPUTE.BACKEND.PARENTS

# missing data imputation with maximum likelihood predictions.
impute.backend.map = function(fitted, data, n, debug = FALSE) {

  nodes = names(data)

  # if there is no missing value, nothing to do.
  missing = !complete.cases(data)

  # call predict() backends so that arguments are not sanitized again.
  for (j in which(missing)) {

    from = nodes[which(!is.na(data[j, ]))]
    to = setdiff(nodes, from)

    # use the obseved part of the observation as the evidence.
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
                  extra = list(from = from, n = n), cluster = NULL, debug = FALSE)

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

    if (debug) {

      cat("  > imputed value:", "\n")
      print(data.frame(estimates), row.names = FALSE)

    }#THEN

    data[j, to] = estimates

  }#FOR

  return(data)

}#IMPUTE.BACKEND.MAP

