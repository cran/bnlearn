
# create a bn.fit object for user-specified local distributions.
custom.fit.backend = function(x, dist, discrete, ordinal, debug = FALSE) {

  # cache node labels.
  nodes = names(x$nodes)
  nnodes = length(nodes)

  # create a dummy bn.fit object from the bn object.
  fitted = structure(vector(nnodes, mode = "list"), names = nodes)
  for (node in nodes)
    fitted[[node]] = list(node = node, parents = x$nodes[[node]]$parents,
                       children = x$nodes[[node]]$children)

  for (node in topological.ordering(x)) {

    # extract the labels of the parents of the node.
    node.parents = x$nodes[[node]]$parents

    if (debug) {

      cat("* processing node", node, ".\n")
      if (length(node.parents) > 0)
        cat("  > found parents:", node.parents, ".\n")

    }#THEN

    # auto-detect the node class from the shape of its parameters:
    if (is.ndmatrix(dist[[node]]))
      param.type = "conditional.probability.table"
    else if (is(dist[[node]], c("lm", "glm", "penfit")))
      param.type = "gaussian.regression"
    else if (is.matrix(dist[[node]]$coef))
      param.type = "conditional.gaussian.regression"
    else if (all(c("inflation", "intensity", "dispersion") %in% names(dist[[node]])))
      param.type = "zero.inflated.hyper.poisson"
    else if (all(c("inflation", "prsucc", "failures") %in% names(dist[[node]])))
      param.type = "zero.inflated.negative.binomial"
    else
      param.type = "gaussian.regression"

    if (param.type == "conditional.probability.table") {

      # sanity check the components of the assignment.
      dist[[node]] = check.dnode.rvalue(dist[[node]], node = node)
      # check the distribution against that of the parents.
      dist[[node]] = check.dnode.rvalue.vs.parents(node, new = dist[[node]],
                       parents = fitted[node.parents])
      # store the new CPT in the bn.fit object.
      fitted[[node]]$prob = normalize.cpt(dist[[node]])
      # set the correct class for method dispatch.
      class(fitted[[node]]) =
        ifelse(node %in% ordinal, "bn.fit.onode", "bn.fit.dnode")

    }#THEN
    else if (param.type == "conditional.gaussian.regression") {

      if (node %in% ordinal)
        stop("node ", node, " is set to be ordinal but is not discrete.")

      # sanity check the components of the assignment.
      dist[[node]] = check.cgnode.rvalue(dist[[node]], node = node)
      # check the assignment against the network structure.
      dist[[node]] =
        check.cgnode.rvalue.vs.parents(node, new = dist[[node]],
          parents = fitted[node.parents])

      # identify discrete and continuous parents and configurations.
      configs = as.character(seq(from = 0, to = ncol(dist[[node]]$coef) - 1))
      # identify discrete and continuous parents.
      fitted[[node]]$dparents = dist[[node]]$dparents
      fitted[[node]]$gparents = dist[[node]]$gparents
      # include the levels of the discrete parents.
      fitted[[node]]$dlevels = dist[[node]]$dlevels
      # store the new coefficients and standard deviations.
      fitted[[node]]$coefficients = noattr(dist[[node]]$coef)
      fitted[[node]]$sd = noattr(dist[[node]]$sd, ok = "names")

      # set the correct class for method dispatch.
      class(fitted[[node]]) = "bn.fit.cgnode"

    }#THEN
    else if (param.type == "gaussian.regression") {

      if (node %in% ordinal)
        stop("node ", node, " is set to be ordinal but is not discrete.")

      # transparently convert regression model objects.
      if (is(dist[[node]], c("lm", "glm", "penfit"))) {

        # ordinary least squares, ridge, lasso, and elastic net.
        dist[[node]] =
          list(coef = .coefficients(dist[[node]]),
               resid = .residuals(dist[[node]]),
               fitted = .fitted(dist[[node]]),
               sd = cgsd(.residuals(dist[[node]]),
                      p = length(.coefficients(dist[[node]]))))

      }#THEN

      # sanity check the components of the assignment.
      dist[[node]] = check.gnode.rvalue(dist[[node]], node = node)
      # check the assignment against the network structure.
      dist[[node]] = check.gnode.rvalue.vs.parents(node, new = dist[[node]],
                       parents = fitted[node.parents])
      # store the new coefficients and standard deviations.
      fitted[[node]]$coefficients = noattr(dist[[node]]$coef, ok = "names")
      fitted[[node]]$sd = noattr(dist[[node]]$sd)
      # set the correct class for method dispatch.
      class(fitted[[node]]) = "bn.fit.gnode"

    }#THEN
    else if (param.type == "zero.inflated.hyper.poisson") {

      if (node %in% ordinal)
        stop("node ", node, " is set to be ordinal but is not discrete.")

      # sanity check the components of the assignment.
      dist[[node]] = check.zihpnode.rvalue(dist[[node]], node = node)
      # check the assignment against the network structure.
      dist[[node]] = check.zihpnode.rvalue.vs.parents(node, new = dist[[node]],
                       parents = fitted[node.parents])
      # store the new coefficients and standard deviations.
      fitted[[node]]$inflation = noattr(dist[[node]]$inflation, ok = "names")
      fitted[[node]]$intensity = noattr(dist[[node]]$intensity, ok = "names")
      fitted[[node]]$dispersion = noattr(dist[[node]]$dispersion)
      # set the correct class for method dispatch.
      class(fitted[[node]]) = "bn.fit.zihpnode"

    }#THEN
    else if (param.type == "zero.inflated.negative.binomial") {

      if (node %in% ordinal)
        stop("node ", node, " is set to be ordinal but is not discrete.")

      # sanity check the components of the assignment.
      dist[[node]] = check.zinbnode.rvalue(dist[[node]], node = node)
      # check the assignment against the network structure.
      dist[[node]] = check.zinbnode.rvalue.vs.parents(node, new = dist[[node]],
                       parents = fitted[node.parents])
      # store the new coefficients and standard deviations.
      fitted[[node]]$inflation = noattr(dist[[node]]$inflation, ok = "names")
      fitted[[node]]$prsucc = noattr(dist[[node]]$prsucc, ok = "names")
      fitted[[node]]$failures = noattr(dist[[node]]$failures)
      # set the correct class for method dispatch.
      class(fitted[[node]]) = "bn.fit.zinbnode"

    }#THEN

   if (debug)
     cat("  > the node has class", class(fitted[[node]]), ".\n")

  }#FOR

  # include fitted values, residuals and configurations if appropriate.
  fitted = bn.fitted.data.values(fitted, dist)
  # guess the correct secondary class ("bn.fit.*net") and return.
  return(structure(fitted, class = c("bn.fit", determine.fitted.class(fitted))))

}#CUSTOM.FIT.BACKEND

# check whether to include fitted values, residuals and configurations.
bn.fitted.data.values = function(fitted, dist) {

  nodes.class = sapply(fitted, class)
  gnodes = names(which(nodes.class == "bn.fit.gnode"))
  cgnodes = names(which(nodes.class == "bn.fit.cgnode"))
  zihpnodes = names(which(nodes.class == "bn.fit.zihpnode"))
  zinbnodes = names(which(nodes.class == "bn.fit.zinbnode"))
  has.values = c(gnodes, cgnodes, zihpnodes, zinbnodes)

  if (length(has.values) > 0) {

    # check whether there is a coherent set of fitted values and residuals.
    nresid = unique(sapply(dist[has.values], function(n) length(n$resid)))
    nfitted = unique(sapply(dist[has.values], function(n) length(n$fitted)))
    nconfig = unique(sapply(dist[cgnodes], function(n) length(n$configs)))

    # all nodes must have residuals and fitted values of the same length.
    full.spec = (length(nresid) == 1) && (length(nfitted) == 1) &&
                    all((nresid > 0) && (nfitted > 0) && (nresid == nfitted))
    # further check discrete parents' configurations for bn.fit.cgnet.
    if (full.spec && (length(cgnodes) > 0))
      full.spec = full.spec && (length(nconfig) == 1) &&
                    all((nconfig > 0) && (nconfig == nfitted))

    # do not trigger a warning if no residuals or fitted values are specified.
    if (!full.spec && (any(nresid > 0) || any(nfitted > 0)))
      warning("different nodes have different number of residuals or fitted values, disregarding.")

  }#THEN

  # cross-check distributions for consistency and populate the bn.fit object.
  for (node in has.values) {

    if (node %in% cgnodes) {

      # save the configurations of the discrete parents.
      if (full.spec)
        fitted[[node]]$configs = noattr(dist[[node]]$configs)
      else
        fitted[[node]]$configs = factor(NA, levels = seq(from = 0,
          to = prod(sapply(fitted[[node]]$dlevels, length)) - 1L))

    }#THEN

    if (full.spec) {

      fitted[[node]]$residuals =
        noattr(dist[[node]]$resid, ok = character(0))
      fitted[[node]]$fitted.values =
        noattr(dist[[node]]$fitted, ok = character(0))

    }#THEN
    else {

      fitted[[node]]$residuals = as.numeric(NA)
      fitted[[node]]$fitted.values = as.numeric(NA)

    }#ELSE

  }#FOR

  return(fitted)

}#BN.FITTED.DATA.VALUES

# determining the network class from the node classes.
determine.fitted.class = function(fitted) {

  nnodes = length(fitted)

  # retrieve and categorize node classes.
  node.classes = factor(sapply(fitted, class), levels = fitted.node.types)
  # count occurrences of each class.
  nct = table(node.classes)
  # scan through valid node class combinations.
  if (nct["bn.fit.dnode"] == nnodes)
    return("bn.fit.dnet")
  else if (nct["bn.fit.onode"] == nnodes)
    return("bn.fit.onet")
  else if (sum(nct[c("bn.fit.dnode", "bn.fit.onode")]) == nnodes)
    return("bn.fit.donet")
  else if (nct["bn.fit.gnode"] == nnodes)
    return("bn.fit.gnet")
  else if (sum(nct[c("bn.fit.gnode", "bn.fit.dnode", "bn.fit.onode",
                     "bn.fit.cgnode")]) == nnodes)
    return("bn.fit.cgnet")
  else if ((nct["bn.fit.zihpnode"] + nct["bn.fit.zinbnode"]) == nnodes)
    return("bn.fit.zinet")
  else
    stop("unsupported combination of node types.")

}#DETERMINE.FITTED.CLASS

