
# estimate the predicted values for a particular node.
predict.bn.fit = function(object, node, data, cluster, method = "parents", ...,
    prob = FALSE, debug = FALSE) {

  # check the data are there.
  data = check.data(data, allow.missing = TRUE, allow.levels = TRUE)
  # a valid node is needed.
  check.nodes(nodes = node, graph = object, max.nodes = 1)
  # check the prediction method.
  check.prediction.method(method, data)
  # check debug and prob.
  check.logical(debug)
  check.logical(prob)

  # check the cluster.
  cluster = check.cluster(cluster)

  if (!is.null(cluster)) {

    # set up the slave processes.
    slaves.setup(cluster)
    # disable debugging, the slaves do not cat() here.
    if (debug) {

      warning("disabling debugging output for parallel computing.")
      debug = FALSE

    }#THEN

  }#THEN

  if (prob && !is(object, c("bn.fit.dnet", "bn.fit.onet", "bn.fit.donet")))
    stop("prediction probabilities are only available for discrete networks.")

  # check optional arguments.
  extra.args = check.prediction.extra.args(method, list(...), node = node,
                 fitted = object, data = data)

  if (method == "parents") {

    # check the fitted model (parents are the only nodes that are actually
    # needed, but checking the target node itself is good if its present).
    if (node %in% names(data))
      nodes.to.check = c(node, object[[node]]$parents)
    else
      nodes.to.check = object[[node]]$parents

  }#THEN
  else if (method == "bayes-lw") {

    # check the fitted model, the conditioning variables and the target node if
    # it is present in the data.
    if (node %in% names(data))
      nodes.to.check = c(node, extra.args$from)
    else
      nodes.to.check = extra.args$from

  }#THEN
  else if (method == "exact") {

    # check the fitted model, the conditioning variables and the target node if
    # it is present in the data.
    if (node %in% names(data))
      nodes.to.check = c(node, extra.args$from)
    else
      nodes.to.check = extra.args$from

    # try to reduce the number of predictors, for speed, but only if all
    # predictors are complete (imputation code takes care of that for incomplete
    # data).
    complete.nodes = attr(data, "metadata")$complete.nodes
    complete.predictors = complete.nodes[extra.args$from]

    if (all(complete.predictors))
      extra.args$from = reduce.predictors.for.exact.inference(fitted = object,
                          target = node, predictors = extra.args$from)

  }#THEN

  # check that the data and the network are consistent.
  check.fit.vs.data(fitted = object, data = data, subset = nodes.to.check)

  predict.backend(fitted = object, node = node, data = data, cluster = cluster,
    method = method, extra.args = extra.args, prob = prob, debug = debug)

}#PREDICT.BN.FIT

# estimate the predicted values for a naive Bayes classfier.
predict.bn.naive = function(object, data, prior, ..., prob = FALSE, debug = FALSE) {

  # check the data are there.
  data = check.data(data, allowed.types = discrete.data.types,
           allow.levels = TRUE)
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

  # get the response variable.
  training = attr(fitted, "training")
  # check the fitted model.
  check.fit.vs.data(fitted = fitted, data = data,
    subset = setdiff(names(fitted), training))
  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  # check the prior distribution.
  prior = check.classifier.prior(prior, fitted[[training]])

  # compute the predicted values.
  naive.classifier(training = training, fitted = fitted, data = data,
    prior = prior, prob = prob, debug = debug)

}#PREDICT.BN.NAIVE

# estimate the predicted values for a TAN classfier.
predict.bn.tan = predict.bn.naive

