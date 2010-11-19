
# passthrough for bn objects.
predict.bn  = function(object, node, data, ...) {

  predict.bn.fit(object = bn.fit(object, data), node = node, data = data, ...)

}#PREDICT.BN

# estimate the predicted values for a particular node.
predict.bn.fit = function(object, node, data, ...) {

  # check the data are there.
  check.data(data)
  # a valid node is needed.
  check.nodes(nodes = node, graph = object, max.nodes = 1)
  # check the fitted model.
  check.fit.vs.data(fitted = object, data = data)
  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  if (is.fitted.discrete(object))
    discrete.prediction(node = node, fitted = object, data = data)
  else
    gaussian.prediction(node = node, fitted = object, data = data)

}#PREDICT.BN.FIT

# estimate the predicted values for a gaussian node.
predict.bn.fit.gnode = function(object, data, ...) {

  nodes = names(data)
  target = object$node

  # check the data are there.
  check.data(data)
  # a valid node is needed.
  check.nodes(nodes = target, graph = nodes, max.nodes = 1)
  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  # create a dummy bn.fit object to pass to gaussian.prediction().
  dummy = vector(length(nodes), mode = "list")
  names(dummy) = nodes
  dummy[[target]] = object
  # compute the predicted values.
  gaussian.prediction(node = target, fitted = dummy, data = data)

}#PREDICT.BN.FIT.GNODE

# estimate the predicted values for a discrete node.
predict.bn.fit.dnode = function(object, data, ...) {

  # check the data are there.
  check.data(data)
  # a valid node is needed.
  check.nodes(nodes = target, graph = nodes, max.nodes = 1)
  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  nodes = names(data)
  target = object$node

  # create a dummy bn.fit object to pass to discrete.prediction().
  dummy = vector(length(nodes), mode = "list")
  names(dummy) = nodes
  dummy[[target]] = object
  # compute the predicted values.
  discrete.prediction(node = target, fitted = dummy, data = data)

}#PREDICT.BN.FIT.DNODE

# estimate the predicted values for a naive Bayes classfier.
predict.bn.naive = function(object, data, prior, ...) {

  # check the data are there.
  check.data(data)
  # check the bn.naive object.
  check.bn.naive(object)

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
  training = root.leaf.nodes(fitted, leaf = FALSE)
  # check the prior distribution.
  prior = check.prior(prior, data[, training])
  # compute the predicted values.
  naive.classifier(training = training, fitted = fitted, data = data,
    prior = prior)

}#PREDICT.BN.NAIVE
