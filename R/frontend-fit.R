
# fit the parameters of the bayesian network for a given network stucture.
bn.fit = function(x, data, method = "mle", ..., debug = FALSE) {

  # check x's class.
  check.bn(x)
  # check the data.
  check.data(data)
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

# replace one conditional probability distribution in a bn.fit object.
"[[<-.bn.fit" = function(x, name, value) {

  # check the label of the node to replace.
  check.nodes(name, x)
  # preserve the original object for subsequent sanity checks.
  to.replace = x[[name]]
  new = to.replace

  if (is(to.replace, "bn.fit.dnode")) {

    # check the consistency of the new conditional distribution.
    value = check.fit.dnode.spec(value)
    # sanity check the new obejct by comparing it to the old one.
    value = check.dnode.vs.spec(value, to.replace)
    # replace the conditional probability table.
    new$prob = value

  }#THEN
  else if (is(to.replace, "bn.fit.gnode")) {

    if (is(value, c("lm", "glm", "penfit"))) {

      # ordinary least squares, ridge, lasso, and elastic net.
      value = list(coef = coefficients(value), resid = residuals(value),
                fitted = fitted(value), sd = sd(residuals(value)))

    }#THEN
    else {

      # check the consistency of the new conditional distribution.
      check.fit.gnode.spec(value)

    }#ELSE

    # sanity check the new obejct by comparing it to the old one.
    check.gnode.vs.spec(value, to.replace)

    # replace the regression coefficients, keeping the names and the ordering.
    if (is.null(names(value$coef)))
      new$coefficients = structure(value$coef, names = names(new$coefficients))
    else
      new$coefficients = value$coef[names(new$coefficients)]

    # replace the residuals' standard deviation.
    if (is.null(value$sd))
      new$sd = sd(value$resid)
    else
      new$sd = value$sd

    # replace the residuals, padding with NAs if needed.
    if (is.null(value$resid))
      new$residuals = rep(as.numeric(NA), length(new$resid))
    else
      new$residuals = value$resid

    # replace the fitted values, padding with NAs if needed.
    if (is.null(value$fitted))
      new$fitted.values = rep(as.numeric(NA), length(new$fitted))
    else
      new$fitted.values = value$fitted

  }#ELSE

  x[name] = list(new)

  return(x)

}#[[<-.BN.FIT

# this is for consistency.
"$<-.bn.fit" = function(x, name, value) {

  `[[<-.bn.fit`(x, name, value)

}#$<-.BN.FIT

custom.fit = function(x, dist) {

  # check x's class.
  check.bn(x)

  nodes = names(x$nodes)
  nnodes = length(nodes)

  # do some basic sanity checks on dist.
  if (!is.list(dist) || is.null(names(dist)))
    stop("the conditional probability distributions must be in a names list.")
  if (length(dist) != nnodes)
    stop("wrong number of conditional probability distributions.")
  check.nodes(names(dist), nodes, min.nodes = nnodes)

  # if all the conditional probability distributions are tables (tables, 
  # matrices and multidimensional are fine), it's a discrete BN.
  discrete = all(sapply(dist, is.ndmatrix))

  # create a dummy bn.fit object from the bn one.
  fitted = structure(vector(nnodes, mode = "list"), names = nodes)

  for (node in nodes) {

    fitted[[node]] = list(node = node, parents = x$nodes[[node]]$parents,
                       children = x$nodes[[node]]$children)

  }#FOR

  if (discrete) {

    # check the consistency of the conditional probability distributions.
    for (cpd in names(dist))
      dist[[cpd]] = check.fit.dnode.spec(dist[[cpd]])

    # cross-check the levels of each node across all CPTs.
    cpt.levels = lapply(dist, function(x) dimnames(x)[[1]])

    for (cpd in names(dist)) {

      # sanity check the new object by comparing it to the old one.
      dist[[cpd]] = check.dnode.vs.spec(dist[[cpd]], old = fitted[[cpd]]$parents,
                      node = cpd, cpt.levels = cpt.levels)
      # store the new CPT in the bn.fit object.
      fitted[[cpd]]$prob = dist[[cpd]]
      # set the correct class for methods' dispatch.
      class(fitted[[cpd]]) = "bn.fit.dnode"

    }#FOR

  }#THEN
  else {

    # convert any lm-type object to the basic list format.
    for (cpd in names(dist)) {

      if (is(dist[[cpd]], c("lm", "glm", "penfit"))) {

        # ordinary least squares, ridge, lasso, and elastic net.
        dist[[cpd]] = list(coef = coefficients(dist[[cpd]]),
                           resid = residuals(dist[[cpd]]),
                           fitted = fitted(dist[[cpd]]),
                           sd = sd(residuals(dist[[cpd]])))

      }#THEN

    }#FOR

    # check whether there is a coherent set of fitted values and residuals.
    nresid = unique(sapply(dist, function(x) length(x$resid)))
    nfitted = unique(sapply(dist, function(x) length(x$fitted)))

    if ((length(nresid) != 1) || (length(nfitted) != 1) || any(nresid != nfitted)) {

      full.spec = FALSE
      warning("different nodes have different number of residuals or fitted values, disregarding.")

    }#THEN
    else {

      full.spec = TRUE

    }#ELSE

    for (cpd in names(dist)) {

      # check the consistency of the conditional probability distribution.
      check.fit.gnode.spec(dist[[cpd]])
      # sanity check the new object by comparing it to the old one.
      check.gnode.vs.spec(dist[[cpd]], old = fitted[[cpd]]$parents,
        node = cpd)
      # store the new CPT in the bn.fit object.
      fitted[[cpd]]$coefficients = dist[[cpd]]$coef
      fitted[[cpd]]$sd = dist[[cpd]]$sd

      if (full.spec) {

        fitted[[cpd]]$residuals = dist[[cpd]]$resid
        fitted[[cpd]]$fitted.values = dist[[cpd]]$fitted

      }#THEN
      else {

        fitted[[cpd]]$residuals = as.numeric(NA)
        fitted[[cpd]]$fitted.values = as.numeric(NA)

      }#ELSE

      # set the correct class for methods' dispatch.
      class(fitted[[cpd]]) = "bn.fit.gnode"

    }#FOR

  }#ELSE

  return(structure(fitted, class = "bn.fit"))

}#CUSTOM.FIT

