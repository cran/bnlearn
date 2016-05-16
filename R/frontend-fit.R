
# fit the parameters of the bayesian network for a given network stucture.
bn.fit = function(x, data, method = "mle", ..., debug = FALSE) {

  # check x's class.
  check.bn(x)
  # check the data.
  if (is(x, c("bn.naive", "bn.tan")))
    check.data(data, allowed.types = discrete.data.types)
  else
    check.data(data)
  # check whether the data agree with the bayesian network.
  check.bn.vs.data(x, data)
  # no parameters if the network structure is only partially directed.
  if (is.pdag(x$arcs, names(x$nodes)))
    stop("the graph is only partially directed.")
  # also check that the network is acyclic.
  if (!is.acyclic(x$arcs, names(x$nodes), directed = TRUE))
    stop("the graph contains cycles.")
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

  # extract the arcs from the fitted network.
  net = empty.graph.backend(names(x))
  arcs(net) = fit2arcs(x)
  # re-create the set of illegal arcs.
  if (is(x, "bn.fit.cgnet"))
    net$learning$illegal = list.cg.illegal.arcs(names(x), x)

  return(net)

}#BN.NET

# extract residuals from continuous bayesian networks.
residuals.bn.fit = function(object, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  if (!is(object, c("bn.fit.gnet", "bn.fit.cgnet")))
    stop("residuals are not defined for discrete bayesian networks.")

  lapply(object, "[[", "residuals")

}#RESIDUALS.BN.FIT

# extract residuals from continuous nodes.
residuals.bn.fit.cgnode = residuals.bn.fit.gnode = function(object, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  object$residuals

}#RESIDUALS.BN.FIT.GNODE

# no residuals here, move along ...
residuals.bn.fit.onode = residuals.bn.fit.dnode = function(object, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  stop("residuals are not defined for discrete nodes.")

}#RESIDUALS.BN.FIT.DNODE

# extract standard errors from continuous bayesian networks.
sigma.bn.fit = function(object, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  if (!is(object, c("bn.fit.gnet", "bn.fit.cgnet")))
    stop("standard errors are not defined for discrete bayesian networks.")

  ll = lapply(object, "[[", "sd")

  # in a conditional Gaussian network, set the standard errors of discrete
  # nodes to NA.
  if (is(object, "bn.fit.cgnet"))
    ll[sapply(ll, is.null)] = NA

  return(ll)

}#SIGMA.BN.FIT

# extract standard errors for continuous nodes.
sigma.bn.fit.cgnode = sigma.bn.fit.gnode = function(object, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  object$sd

}#SIGMA.BN.FIT.GNODE

# no sigma here, move along ...
sigma.bn.fit.onode = sigma.bn.fit.dnode = function(object, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  stop("standard errors are not defined for discrete nodes.")

}#SIGMA.BN.FIT.DNODE

# extract fitted values from continuous bayesian networks.
fitted.bn.fit = function(object, ...) {

  if (!is(object, c("bn.fit.gnet", "bn.fit.cgnet")))
    stop("fitted values are not defined for discrete bayesian networks.")

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  lapply(object, "[[", "fitted.values")

}#FITTED.BN.FIT

# extract fitted values from continuous nodes.
fitted.bn.fit.cgnode = fitted.bn.fit.gnode = function(object, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  object$fitted.values

}#FITTED.BN.FIT.GNODE

# no fitted values here, move along ...
fitted.bn.fit.onode = fitted.bn.fit.dnode = function(object, ...) {

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
coef.bn.fit.cgnode = coef.bn.fit.gnode = function(object, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  object$coefficients

}#COEF.BN.FIT.GNODE

# extract probabilities from discrete nodes.
coef.bn.fit.onode = coef.bn.fit.dnode = function(object, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  object$prob

}#COEF.BN.FIT.DNODE

# logLik method for class 'bn.fit'.
logLik.bn.fit = function(object, data, nodes, by.sample = FALSE, ...) {

  # check the data are there.
  check.data(data)
  # check the fitted model.
  check.fit.vs.data(fitted = object, data = data)
  # warn about unused arguments.
  check.unused.args(list(...), character(0))
  # check the nodes whose logLik components we are going to compute.
  if (missing(nodes))
    nodes = names(object)
  else
    check.nodes(nodes, object)

  llik = entropy.loss(fitted = object, data = data, keep = nodes,
           by.sample = by.sample)$loss

  if (!by.sample)
    llik = - nrow(data) * llik

  return(llik)

}#LOGLIK.BN.FIT

# AIC method for class 'bn.fit'.
AIC.bn.fit = function(object, data, ..., k = 1) {

  logLik(object, data) - k * nparams(object)

}#AIC.BN.FIT

# BIC method for class 'bn.fit'.
BIC.bn.fit = function(object, data, ...) {

  logLik(object, data) - log(nrow(data))/2 * nparams(object)

}#BIC.BN.FIT

# replace one conditional probability distribution in a bn.fit object.
"[[<-.bn.fit" = function(x, name, value) {

  # check x's class.
  check.fit(x)
  # check the label of the node to replace.
  check.nodes(name, x)

  x[name] = list(fitted.assignment.backend(x, name, value))

  return(x)

}#[[<-.BN.FIT

# this is for consistency.
"$<-.bn.fit" = function(x, name, value) {

  `[[<-.bn.fit`(x, name, value)

}#$<-.BN.FIT

# create a bn.fit object for user-specified local distributions.
custom.fit = function(x, dist, ordinal) {

  # check x's class.
  check.bn(x)
  # cache node labels.
  nodes = names(x$nodes)
  nnodes = length(nodes)

  # no parameters if the network structure is only partially directed.
  if (is.pdag(x$arcs, nodes))
    stop("the graph is only partially directed.")
  # also check that the network is acyclic.
  if (!is.acyclic(x$arcs, nodes, directed = TRUE))
    stop("the graph contains cycles.")

  # do some basic sanity checks on dist.
  if (!is.list(dist) || is.null(names(dist)))
    stop("the conditional probability distributions must be in a names list.")
  if (length(dist) != nnodes)
    stop("wrong number of conditional probability distributions.")
  check.nodes(names(dist), nodes, min.nodes = nnodes)

  # only discrete nodes can be parameterized by a CPT, all the others are lists
  # with multiple elements.
  discrete = sapply(dist, is.ndmatrix)

  # check ordinal ...
  if (missing(ordinal))
    ordinal = character(0)
  else
    check.nodes(ordinal, graph = nodes)
  # ... and that nodes that are supposed to be ordinal are discrete in the
  # first place.
  if (!all(discrete[ordinal]))
    stop("node(s)", paste0(" '", names(discrete[!discrete[ordinal]]), "'"),
      " are set to be ordinal but are not discrete.")

  custom.fit.backend(x = x, dist = dist, discrete = discrete, ordinal = ordinal)

}#CUSTOM.FIT

