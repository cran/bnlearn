
# fit the parameters of the bayesian network for a given network stucture.
bn.fit = function(x, data, method = "mle", ..., debug = FALSE) {

  # check x's class.
  check.bn(x)
  # check whether the data agree with the bayesian network.
  check.bn.vs.data(x, data)
  # no parameters if the network structure is only partially directed.
  if (is.pdag(x$arcs, names(x$nodes)))
    stop("the graph is only partially directed.")
  # check the fitting method (maximum likelihood, bayesian, etc.)
  check.fitting.method(method, data)
  # check the extra arguments.
  extra.args = check.fitting.args(method, x, data, list(...))

  bn.fit.backend(x = x, data = data, method = method, extra.args = extra.args,
    debug = debug)

}#BN.FIT

# get back the network structure from the fitted object.
bn.net = function(x, debug = FALSE) {

  # check x's class.
  check.fit(x)

  nodes = names(x)
  net = empty.graph.backend(nodes)
  arcs(net) = fit2arcs(x)

  return(net)

}#BN.NET

# extract residuals for continuous bayesian networks.
residuals.bn.fit = function(object, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  if (is.fitted.discrete(object))
    stop("residuals are not defined for discrete bayesian networks.")

  lapply(object, "[[", "residuals")

}#RESIDUALS.BN.FIT

# extract residuals for continuous nodes.
residuals.bn.fit.gnode = function(object, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  object$residuals

}#RESIDUALS.BN.FIT.GNODE

# no residuals here, move along ...
residuals.bn.fit.dnode = function(object, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  stop("residuals are not defined for discrete nodes.")

}#RESIDUALS.BN.FIT.DNODE

# extract fitted values for continuous bayesian networks.
fitted.bn.fit = function(object, ...) {

  if (is.fitted.discrete(object))
    stop("fitted values are not defined for discrete bayesian networks.")

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  lapply(object, "[[", "fitted.values")

}#FITTED.BN.FIT

# extract fitted values for continuous nodes.
fitted.bn.fit.gnode = function(object, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  object$fitted.values

}#FITTED.BN.FIT.GNODE

# no fitted values here, move along ...
fitted.bn.fit.dnode = function(object, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  stop("fitted values are not defined for discrete nodes.")

}#FITTED.BN.FIT.DNODE

# extract parameters from any bayesian network.
coef.bn.fit = function(object, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  lapply(object, "coef")

}#COEF.BN.FIT

# extract regression coefficients from continuous nodes.
coef.bn.fit.gnode = function(object, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  object$coefficients

}#COEF.BN.FIT.GNODE

# extract probabilities from discrete nodes.
coef.bn.fit.dnode = function(object, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  object$prob

}#COEF.BN.FIT.DNODE

# logLik method for class 'bn.fit'.
logLik.bn.fit = function(object, data, ...) {

  # check the data are there.
  check.data(data)
  # check the fitted model.
  check.fit.vs.data(fitted = object, data = data)
  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  nodes = names(object)
  ndata = nrow(data)

  # parameter sanitization done in the score() function.
  if (is.data.discrete(data))
    - ndata * discrete.loss(nodes = nodes, fitted = object, data = data)$loss
  else
    - ndata * gaussian.loss(nodes = nodes, fitted = object, data = data)$loss

}#LOGLIK.BN.FIT

# AIC method for class 'bn.fit'.
AIC.bn.fit = function(object, data, ..., k = 1) {

  logLik(object, data) - k * nparams(object)

}#AIC.BN.FIT

