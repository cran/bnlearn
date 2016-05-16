
# passthrough for bn objects.
predict.bn  = function(object, node, data, method = "parents", ...,
    prob = FALSE, debug = FALSE) {

  predict.bn.fit(object = bn.fit(object, data), node = node, data = data,
    method = method, ..., prob = prob, debug = debug)

}#PREDICT.BN

# estimate the predicted values for a particular node.
predict.bn.fit = function(object, node, data, method = "parents", ...,
    prob = FALSE, debug = FALSE) {

  # check the data are there.
  check.data(data, allow.levels = TRUE)
  # a valid node is needed.
  check.nodes(nodes = node, graph = object, max.nodes = 1)
  # check the prediction method.
  check.prediction.method(method, data)
  # check debug and prob.
  check.logical(debug)
  check.logical(prob)

  if (prob && !is(object, c("bn.fit.dnet", "bn.fit.onet", "bn.fit.donet")))
    stop("prediction probabilities are only available for discrete networks.")

  # warn about unused arguments.
  extra.args = list(...)
  check.unused.args(extra.args, prediction.extra.args[[method]])

  if (method == "parents") {

    # check the fitted model (parents are the only nodes that are actually
    # needed).
    check.fit.vs.data(fitted = object, data = data,
      subset = object[[node]]$parents)

    if (is(object, c("bn.fit.dnet", "bn.fit.onet", "bn.fit.donet")))
      discrete.prediction(node = node, fitted = object, data = data,
        prob = prob, debug = debug)
    else if (is(object, "bn.fit.gnet"))
      gaussian.prediction(node = node, fitted = object, data = data,
        debug = debug)
    else if (is(object, "bn.fit.cgnet"))
      mixedcg.prediction(node = node, fitted = object, data = data,
        debug = debug)

  }#THEN
  else if (method == "bayes-lw") {

    # check the variables to predict from.
    if (is.null(extra.args$from))
      extra.args$from = setdiff(names(data), node)
    else
      check.nodes(nodes = extra.args$from, graph = object, min.nodes = 1)
    # check that they do not include the node to predict.
    if (node %in% extra.args$from)
      stop("node ", node, " is both a predictor and being predicted.")
    # check the fitted model and the conditioning variables.
    check.fit.vs.data(fitted = object, data = data, subset = extra.args$from)
    # check the number of particles to be used for each prediction.
    if (is.null(extra.args$n))
      extra.args$n = 500
    else if (!is.positive.integer(extra.args$n))
      stop("the number of observations to be sampled must be a positive integer number.")

    map.prediction(node = node, fitted = object, data = data, n = extra.args$n,
      from = extra.args$from, prob = prob, debug = debug)

  }#THEN

}#PREDICT.BN.FIT

# estimate the predicted values for a naive Bayes classfier.
predict.bn.naive = function(object, data, prior, ..., prob = FALSE, debug = FALSE) {

  # check the data are there.
  check.data(data, allowed.types = discrete.data.types, allow.levels = TRUE)
  # check the bn.{naive,tan} object.
  if (is(object, "bn.naive"))
    check.bn.naive(object)
  else
    check.bn.tan(object)

  # check debug and prob.
  check.logical(debug)
  check.logical(prob)

  # fit the network if needed.
  if (is(object, "bn"))
    fitted = bn.fit(object, data)
  else
    fitted = object

  # check the fitted model.
  check.fit.vs.data(fitted = fitted, data = data)
  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  # get the response variable.
  training = attr(fitted, "training")
  # check the prior distribution.
  prior = check.classifier.prior(prior, fitted[[training]])

  # compute the predicted values.
  naive.classifier(training = training, fitted = fitted, data = data,
    prior = prior, prob = prob, debug = debug)

}#PREDICT.BN.NAIVE

# estimate the predicted values for a TAN classfier.
predict.bn.tan = predict.bn.naive

