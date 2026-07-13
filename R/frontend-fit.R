
# fit the parameters of the bayesian network for a given network structure.
bn.fit = function(x, data, cluster, method, ..., keep.fitted = TRUE,
    debug = FALSE) {

  check.bn(x)
  # check the data.
  if (is(x, available.classifiers)) {

    data = check.data(data, allow.missing = TRUE,
             allowed.types = discrete.data.types)

  }#THEN
  else {

    data = check.data(data, allow.missing = TRUE)

  }#ELSE
  # check the fitting method (maximum likelihood, bayesian, etc.)
  method = check.fitting.method(method, data)
  # check whether the data agree with the bayesian network.
  data = check.bn.vs.data(x, data, reorder = grepl("hard-em", method))
  # no parameters if the network structure is only partially directed.
  if (!is.completely.directed(x))
    stop("the graph is only partially directed.")
  # also check that the network is acyclic.
  if (!is.acyclic(x$arcs, names(x$nodes), directed = TRUE))
    stop("the graph contains cycles.")
  # check whether the network is valid for the method.
  check.arcs.against.assumptions(x$arcs, data, method)
  # check the extra arguments.
  extra.args = check.fitting.args(method, x, data, list(...))

  check.logical(debug)
  check.logical(keep.fitted)

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

  bn.fit.backend(x = x, data = data, cluster = cluster, method = method,
    extra.args = extra.args, keep.fitted = keep.fitted, debug = debug)

}#BN.FIT

# get back the network structure from the fitted object.
bn.net = function(x) {

  check.bn.or.fit(x)

  # nothing to do, the input is already a network structure().
  if (is(x, "bn"))
    return(x)

  # extract the arcs from the fitted network.
  net = empty.graph.backend(names(x))
  arcs(net) = fit2arcs(x)
  # re-create the set of illegal arcs.
  type = class(x)
  net$learning$illegal = list.illegal.arcs(names(x), x, type[length(type)])
  # preserve node causal roles, if present.
  attrs = attributes(x)
  if ("roles" %in% names(attrs))
    net$learning$roles = attrs$roles
  # also preserve the training node label set by classifiers.
  if ("training" %in% names(attrs))
    net$learning$args$training = attrs$training
  # also preserve auxiliary classes related to causal inference and classifiers.
  all.classes = setdiff(attrs$class, available.fitted)
  all.classes[all.classes == "bn.fit"] = "bn"
  class(net) = all.classes

  return(net)

}#BN.NET

# extract residuals from continuous bayesian networks.
residuals.bn.fit = function(object, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  if (is(object, c("bn.fit.dnet", "bn.fit.onet", "bn.fit.donet")))
    stop("residuals are not defined for discrete bayesian networks.")

  lapply(object, "[[", "residuals")

}#RESIDUALS.BN.FIT

# extract residuals from continuous and zero-inflated nodes.
residuals.bn.fit.cgnode = residuals.bn.fit.gnode = residuals.bn.fit.zihpnode =
  residuals.bn.fit.zinbnode = function(object, ...) {

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

  if (is(object, c("bn.fit.dnet", "bn.fit.onet", "bn.fit.donet")))
    stop("standard errors are not defined for discrete bayesian networks.")

  ll = lapply(object, "[[", "sd")

  # in a conditional Gaussian network, set the standard errors of discrete
  # nodes to NA.
  if (is(object, "bn.fit.cgnet"))
    ll[sapply(ll, is.null)] = NA

  return(ll)

}#SIGMA.BN.FIT

# extract standard errors for continuous nodes.
sigma.bn.fit.gnode = function(object, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  object$sd

}#SIGMA.BN.FIT.GNODE

sigma.bn.fit.cgnode = function(object, for.parents, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), "for.parents")

  if (missing(for.parents)) {

    sd = object$sd

  }#THEN
  else {

    for.parents = check.discrete.parents.configuration(for.parents, object)

     # enumerate all possible configurations...
    all.configurations = expand.grid(object$dlevels, stringsAsFactors = FALSE)
    # ... find which one to return...
    requested = which(apply(all.configurations, 1, identical, unlist(for.parents)))
    # ... and extract it.
    sd = noattr(object$sd[requested])

  }#ELSE

  return(sd)

}#SIGMA.BN.FIT.CGNODE

# no sigma here, move along ...
sigma.bn.fit.onode = sigma.bn.fit.dnode = function(object, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  stop("standard errors are not defined for discrete nodes.")

}#SIGMA.BN.FIT.DNODE

# no sigma here, move along ...
sigma.bn.fit.zihpnode = sigma.bn.fit.zinbnode = function(object, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  stop("standard errors are not defined for zero-inflated nodes.")

}#SIGMA.BN.FIT.ZIHPNODE

# extract fitted values from continuous bayesian networks.
fitted.bn.fit = function(object, ...) {

  if (is(object, c("bn.fit.dnet", "bn.fit.onet", "bn.fit.donet")))
    stop("fitted values are not defined for discrete bayesian networks.")

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  lapply(object, "[[", "fitted.values")

}#FITTED.BN.FIT

# extract fitted values from continuous and zero-inflated nodes.
fitted.bn.fit.cgnode = fitted.bn.fit.gnode = fitted.bn.fit.zihpnode =
  fitted.bn.fit.zinbnode = function(object, ...) {

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
coef.bn.fit.gnode = function(object, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  object$coefficients

}#COEF.BN.FIT.GNODE

coef.bn.fit.cgnode = function(object, for.parents, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), "for.parents")

  if (missing(for.parents)) {

    coefficients = object$coefficients

  }#THEN
  else {

    for.parents = check.discrete.parents.configuration(for.parents, object)

    # enumerate all possible configurations...
    all.configurations = expand.grid(object$dlevels, stringsAsFactors = FALSE)
    # ... find which one to return...
    requested = which(apply(all.configurations, 1, identical, unlist(for.parents)))
    # ... and extract it.
    coefficients = object$coefficients[, requested]

  }#ELSE

  return(coefficients)

}#COEF.BN.FIT.CGNODE

# extract probabilities from discrete nodes.
coef.bn.fit.onode = coef.bn.fit.dnode = function(object, for.parents, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), "for.parents")

  if (missing(for.parents)) {

    coefficients = object$prob

  }#THEN
  else {

    for.parents = check.discrete.parents.configuration(for.parents, object)

    # reorder the dimensions...
    requested = c(list(TRUE), for.parents[object$parents])
    # ... and extract the probabilities.
    coefficients = do.call("[", c(list(object$prob), requested))

  }#ELSE

  return(coefficients)

}#COEF.BN.FIT.DNODE

# extract the parameters of the zero-inflated nodes.
coef.bn.fit.zihpnode = coef.bn.fit.zinbnode = function(object, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  if (is(object, "bn.fit.zihpnode"))
    object[c("inflation", "intensity", "dispersion")]
  else
    object[c("inflation", "prsucc", "failures")]

}#COEF.BN.FIT.ZIHPNODE

# logLik method for class 'bn.fit'.
logLik.bn.fit = function(object, data, nodes, by.sample = FALSE,
    na.rm = FALSE, debug = FALSE, ...) {

  # check the data are there.
  data = check.data(data, allow.missing = TRUE, allow.levels = TRUE)
  # check the fitted model.
  data = check.fit.vs.data(fitted = object, data = data, reorder = TRUE)
  # warn about unused arguments.
  check.unused.args(list(...), character(0))
  # check the nodes whose logLik components we are going to compute.
  if (missing(nodes))
    nodes = names(object)
  else
    check.nodes(nodes, object)
  # check the logical arguments.
  check.logical(by.sample)
  check.logical(na.rm)
  check.logical(debug)

  loglikelihood(fitted = object, data = data, by.sample = by.sample,
    keep = nodes, propagate.missing = !na.rm, debug = debug)

}#LOGLIK.BN.FIT

# AIC method for class 'bn.fit'.
AIC.bn.fit = function(object, data, ..., k = 1) {

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  ll = logLik.bn.fit(object, data)

  return(ll - k * attr(ll, "df"))

}#AIC.BN.FIT

# BIC method for class 'bn.fit'.
BIC.bn.fit = function(object, data, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  ll = logLik.bn.fit(object, data)

  return(ll - log(nrow(data)) / 2 * attr(ll, "df"))

}#BIC.BN.FIT

# replace one conditional probability distribution in a bn.fit object.
"[[<-.bn.fit" = function(x, name, value) {

  check.fit(x)
  # check the label of the node to replace.
  check.nodes(name, x)

  x[name] = list(bn.fit.assignment.backend(x, name, value))

  return(x)

}#[[<-.BN.FIT

# this is for consistency.
"$<-.bn.fit" = function(x, name, value) {

  `[[<-.bn.fit`(x, name, value)

}#$<-.BN.FIT

# create a bn.fit object for user-specified local distributions.
custom.fit = function(x, dist, ordinal, debug = FALSE) {

  check.bn(x)
  check.logical(debug)
  # cache the node labels.
  nodes = names(x$nodes)
  nnodes = length(nodes)

  # no parameters if the network structure is only partially directed.
  if (!is.completely.directed(x))
    stop("the graph is only partially directed.")
  # also check that the network is acyclic.
  if (!is.acyclic(x$arcs, nodes, directed = TRUE))
    stop("the graph contains cycles.")

  # do some basic sanity checks on dist.
  if (!is.list(dist) || is.null(names(dist)))
    stop("the conditional probability distributions must be in a named list.")
  if (length(dist) != nnodes)
    stop("wrong number of conditional probability distributions.")
  check.nodes(names(dist), nodes, min.nodes = nnodes)

  # check ordinal.
  if (missing(ordinal))
    ordinal = character(0)
  else
    check.nodes(ordinal, graph = nodes)

  custom.fit.backend(x = x, dist = dist, ordinal = ordinal, debug = debug)

}#CUSTOM.FIT

# Shannon's entropy of a fitted network.
H = function(P) {

  # check the bn.fit object.
  check.fit(P)

  if (is(P, "bn.fit.zinet"))
    stop("zero inflated networks are not supported.")

  shannon.entropy(P)

}#H

# Kullback-Leibler divergence between a network Q and a reference network P.
KL = function(P, Q) {

  # both should be fitted networks...
  check.fit(P)
  check.fit(Q)

  if (is(P, "bn.fit.zinet") || is(Q, "bn.fit.zinet"))
    stop("zero inflated networks are not supported.")

  # ... and they should have compatible distributions.
  P = check.fitted.vs.fitted(P, Q, local = FALSE)

  kullback.leibler(P, Q)

}#KL

# does the bn.fit object contain any NA parameter values?
identifiable = function(x, by.node = FALSE) {

  # check the bn.fit object.
  check.fit(x)

  per.node = sapply(x, function(ld) {

    if (is(ld, c("bn.fit.dnode", "bn.fit.onode")))
      return(!anyNA(ld$prob))
    else if (is(ld, c("bn.fit.gnode", "bn.fit.cgnode")))
      !anyNA(ld$coefficients) && !anyNA(ld$sd)
    else if (is(ld, "bn.fit.zihpnode"))
      !anyNA(ld$inflation) && !anyNA(ld$intensity) && !is.na(ld$dispersion)
    else if (is(ld, "bn.fit.zinbnode"))
      !anyNA(ld$inflation) && !anyNA(ld$prsucc) && !is.na(ld$failures)
    else
      stop("unknown node type '", class(ld), "'.")

  })

  # the model is identifiable if all nodes are identifiable.
  if (by.node)
    return(per.node)
  else
    return(all(per.node))

}#IDENTIFIABLE

# does the bn.fit object encode a singular model?
singular = function(x, by.node = FALSE) {

  # check the bn.fit object.
  check.fit(x)
  if (is(x, "bn.fit.zinet"))
    stop("zero inflated networks are not supported.")

  per.node = sapply(x, function(ld) {

    if (is(ld, c("bn.fit.dnode", "bn.fit.onode"))) {

      cpt.dims = dim(ld$prob)

      if (length(cpt.dims) == 1)
        return(all(ld$prob %in% c(0, 1)))
      else {

        # cast the CPT into two dimensions, so that we can ...
        two.d = array(ld$prob, dim = c(cpt.dims[1], prod(cpt.dims[-1])))
        # ... check all conditional distributions in a single apply() call.
        sing = apply(two.d, 2, function(p) all(p %in% c(0, 1)))
        # discount NAs from unidentifiable conditional distributions.
        return(any(sing, na.rm = TRUE))

      }#ELSE

    }#THEN
    else if (is(ld, c("bn.fit.gnode", "bn.fit.cgnode"))) {

      return(any(ld$sd == 0, na.rm = TRUE))

    }#THEN

  })

  # the model is singular if at least one node is singular.
  if (by.node)
    return(per.node)
  else
    return(any(per.node))

}#SINGULAR
