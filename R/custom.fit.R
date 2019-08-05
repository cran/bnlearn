
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

    if (is.ndmatrix(dist[[node]])) {

      # first self-check the local distribution.
      dist[[node]] = check.dnode(dist[[node]], node = node)
      # check thee distribution against that of the parents.
      dist[[node]] = check.dnode.vs.parents(node, new = dist[[node]],
                       parents = fitted[node.parents])
      # store the new CPT in the bn.fit object.
      fitted[[node]]$prob = normalize.cpt(dist[[node]])
      # set the correct class for method dispatch.
      class(fitted[[node]]) =
        ifelse(node %in% ordinal, "bn.fit.onode", "bn.fit.dnode")

    }#THEN
    else {

      if (node %in% ordinal)
        stop("node", node, " is set to be ordinal but are not discrete.")

      # transparently convert regression models' objects.
      if (is(dist[[node]], c("lm", "glm", "penfit"))) {

        # ordinary least squares, ridge, lasso, and elastic net.
        dist[[node]] =
          list(coef = .coefficients(dist[[node]]),
               resid = .residuals(dist[[node]]),
               fitted = .fitted(dist[[node]]),
               sd = cgsd(.residuals(dist[[node]]),
                      p = length(.coefficients(dist[[node]]))))

      }#THEN

      dist[[node]] = check.gnode(dist[[node]], node = node)
      # sanity check the distribution by comparing it to the network structure.
      if (is(dist[[node]]$coef, "matrix")) {

        dist[[node]] =
          check.cgnode.vs.parents(node, new = dist[[node]],
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
      else {

        dist[[node]] = check.gnode.vs.parents(node, new = dist[[node]],
                         parents = fitted[node.parents])
        # store the new coefficients and standard deviations.
        fitted[[node]]$coefficients = noattr(dist[[node]]$coef, ok = "names")
        fitted[[node]]$sd = noattr(dist[[node]]$sd)
        # set the correct class for method dispatch.
        class(fitted[[node]]) = "bn.fit.gnode"

      }#ELSE

    }#ELSE

   if (debug)
     cat("  > the node has class", class(fitted[[node]]), ".\n")

  }#FOR

  # include fitted values, residuals and configurations if appropriate.
  fitted = full.spec.backend(fitted, dist)
  # guess the correct secondary class ("bn.fit.*net") and return.
  return(structure(fitted, class = c("bn.fit", determine.fitted.class(fitted))))

}#CUSTOM.FIT.BACKEND

# check whether to include fitted values, residuals and configurations.
full.spec.backend = function(fitted, dist) {

  nodes.class = sapply(fitted, class)
  gnodes = names(nodes.class[nodes.class == "bn.fit.gnode"])
  cgnodes = names(nodes.class[nodes.class == "bn.fit.cgnode"])
  both = c(gnodes, cgnodes)

  if (length(both) > 0) {

    # check whether there is a coherent set of fitted values and residuals.
    nresid = unique(sapply(dist[both], function(n) length(n$resid)))
    nfitted = unique(sapply(dist[both], function(n) length(n$fitted)))
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
  for (node in both) {

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

}#FULL.SPEC.BACKEND

# return the corect class based on the node classes.
determine.fitted.class = function(fitted) {

  nnodes = length(fitted)

  # retrieve and classify node classes.
  node.classes = factor(sapply(fitted, class), levels = fitted.node.types)
  # count occurrences of each class.
  nct = table(node.classes)

  if (nct["bn.fit.dnode"] == nnodes)
    return("bn.fit.dnet")
  else if (nct["bn.fit.onode"] == nnodes)
    return("bn.fit.onet")
  else if (nct["bn.fit.gnode"] == nnodes)
    return("bn.fit.gnet")
  else if ((nct["bn.fit.gnode"] == 0) && (nct["bn.fit.cgnode"] == 0))
    return("bn.fit.donet")
  else
    return("bn.fit.cgnet")

}#GUESS.FITTED.CLASS

