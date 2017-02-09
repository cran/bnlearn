
# missing data imputation with maximum likelihood predictions.
impute.backend.parents = function(fitted, data, debug = FALSE) {

  # check the variables in topological order, to ensure parents are complete.
  for (i in schedule(fitted)) {

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

  nodes = nodes(fitted)

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
      evidence = lapply(data[j, from], 
                   function(x) if (is.factor(x)) as.character(x))

    # simulate the particles and the weights using likelihood weighting.
    particles = 
      conditional.probability.query(fitted = fitted, event = to,
        evidence = evidence, method = "lw",
        extra = list(from = from, n = n), probability = FALSE, 
        cluster = NULL, debug = FALSE)

    # impute by posterior mode (discrete variables) or posterior expectation
    # (continuous variables).
    particles = sapply(particles, function(x, w) {

      if (is.factor(x))
        names(which.max(by(w, INDICES = x, FUN = sum)))
      else if (is.numeric(x))
        weighted.mean(x = x, w = w)

    }, w = attr(particles, "weights"))

    data[j, to] = particles

  }#FOR

  return(data)

}#IMPUTE.BACKEND.MAP

