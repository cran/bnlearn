
# check the method used for prediction.
check.prediction.method = function(method, data) {

  if (!missing(method) && !is.null(method)) {

    check.label(method, choices = available.prediction.methods,
      labels = prediction.labels, argname = "prediction method", see = "predict")

    return(method)

  }#THEN
  else {

    return("parents")

  }#ELSE

}#CHECK.PREDICTION.METHOD

# sanitize the extra arguments passed to prediction methods.
check.prediction.extra.args = function(method, extra.args, node, fitted, data) {

  # check the number of particles to be used for each prediction.
  if (has.argument(method, "n", prediction.extra.args)) {

    if (is.null(extra.args[["n"]]))
      extra.args[["n"]] = 500
    else if (!is.positive.integer(extra.args[["n"]]))
      stop("the number of observations to be sampled must be a positive integer number.")

  }#THEN

  # check labels of the nodes to predict from.
  if (has.argument(method, "from", prediction.extra.args)) {

    # check the variables to predict from, using the network as a reference.
    if (is.null(extra.args[["from"]]))
      extra.args$from = intersect(setdiff(names(data), node), names(fitted))
    else
      check.nodes(nodes = extra.args[["from"]], graph = fitted, min.nodes = 0)
    # check that they do not include the node to predict.
    if (node %in% extra.args[["from"]])
      stop("node ", node, " is both a predictor and being predicted.")

  }#THEN

  # warn about and remove unused arguments.
  extra.args = check.unused.args(extra.args, prediction.extra.args[[method]])

  return(extra.args)

}#CHECK.PREDICTION.EXTRA.ARGS

# check the method used for imputation.
check.imputation.method = function(method, data) {

  if (!missing(method) && !is.null(method)) {

    check.label(method, choices = available.imputation.methods,
      labels = imputation.labels, argname = "imputation method", see = "impute")

    return(method)

  }#THEN
  else {

    return("bayes-lw")

  }#ELSE

}#CHECK.IMPUTATION.METHOD

# sanitize the extra arguments passed to imputation methods.
check.imputation.extra.args = function(method, extra.args) {

  # check the number of particles to be used for each imputation.
  if (has.argument(method, "n", imputation.extra.args)) {

    if (is.null(extra.args[["n"]]))
      extra.args[["n"]] = 500
    else if (!is.positive.integer(extra.args[["n"]]))
      stop("the number of observations to be sampled must be a positive integer number.")

  }#THEN

  # warn about and remove unused arguments.
  extra.args = check.unused.args(extra.args, imputation.extra.args[[method]])

  return(extra.args)

}#CHECK.IMPUTATION.EXTRA.ARGS

# remove redundant predictors that are d-separated from the target node by other
# predictors.
reduce.predictors.for.exact.inference = function(fitted, target, predictors) {

  markov.blanket = mb.fitted(fitted, target)
  if (all(markov.blanket %in% predictors)) {

    # there is no point in predicting from a superset of the Markov blanket of
    # the target node, just use the Markov blanket.
    predictors = markov.blanket

  }#THEN
  else if (length(predictors) >= 2) {

    dag = bn.net(fitted)

    # see if any predictor is d-separated from the target node by the rest and,
    # if so, drop it.
    for (candidate in predictors) {

      if (length(predictors) == 1)
        break

      reduced = setdiff(predictors, candidate)

      if (dseparation(dag, x = target, y = candidate, z = reduced))
        predictors = reduced

    }#FOR

  }#ELSE

  return(predictors)

}#REDUCE.PREDICTORS.FOR.EXACT.INFERENCE

