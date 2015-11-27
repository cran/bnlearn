
# create a bn.fit object for user-specified local distributions.
custom.fit.backend = function(x, dist, discrete, ordinal) {

  # cache node labels.
  nodes = names(x$nodes)
  nnodes = length(nodes)

  # create a dummy bn.fit object from the bn object.
  fitted = structure(vector(nnodes, mode = "list"), names = nodes)
  for (node in nodes)
    fitted[[node]] = list(node = node, parents = x$nodes[[node]]$parents,
                       children = x$nodes[[node]]$children)

  # first self-check each local distribution.
  nconfig = nresid = nfitted = structure(rep(NA, nnodes), names = nodes)

  for (node in nodes) {

    if (discrete[node]) {

      dist[[node]] = check.fit.dnode.spec(dist[[node]], node = node)
      # set the correct class for method dispatch.
      class(fitted[[node]]) =
        ifelse(node %in% ordinal, "bn.fit.onode", "bn.fit.dnode")

    }#THEN
    else {

      # transparently convert regression models' objects.
      if (is(dist[[node]], c("lm", "glm", "penfit"))) {

        # ordinary least squares, ridge, lasso, and elastic net.
        dist[[node]] =
          list(coef = minimal.coefficients(dist[[node]]),
               resid = minimal.residuals(dist[[node]]),
               fitted = minimal.fitted(dist[[node]]),
               sd = cgsd(minimal.residuals(dist[[node]]),
                      p = length(minimal.coefficients(dist[[node]]))))

      }#THEN

      dist[[node]] = check.fit.gnode.spec(dist[[node]], node = node)
      # sanity check the distribution by comparing it to the network structure.
      if (is(dist[[node]]$coef, "matrix")) {

        check.cgnode.vs.spec(dist[[node]], old = fitted[[node]]$parents,
          node = node)
        # set the correct class for method dispatch.
        class(fitted[[node]]) = "bn.fit.cgnode"

        if ("configs" %in% names(dist[[node]]))
          nconfig[node] = length(dist[[node]]$configs)

      }#THEN
      else {

        dist[[node]] = check.gnode.vs.spec(dist[[node]],
                         old = fitted[[node]]$parents, node = node)
        # set the correct class for method dispatch.
        class(fitted[[node]]) = "bn.fit.gnode"

      }#ELSE

      if ("resid" %in% names(dist[[node]]))
        nresid[node] = length(dist[[node]]$resid)
      if ("fitted" %in% names(dist[[node]]))
        nfitted[node] = length(dist[[node]]$fitted)

    }#ELSE

  }#FOR

  # guess the correct secondary class ("bn.fit.*net").
  secondary.class = guess.fitted.class(fitted)

  # check whether there is a coherent set of fitted values and residuals.
  if (!all(discrete)) {

    nresid = unique(nresid[!discrete])
    nfitted = unique(nfitted[!discrete])
    nconfig = unique(nconfig[sapply(fitted, class) == "bn.fit.cgnode"])

    # all nodes must have residuals and fitted values of the same length.
    full.spec = all(!is.na(nresid) & !is.na(nfitted)) &&
                (length(nresid) == 1) && (length(nfitted) == 1)
    if (full.spec)
      full.spec = full.spec && (nresid == nfitted)
    # further check discrete parents' configurations for bn.fit.cgnet.
    if (full.spec && (secondary.class == "bn.fit.cgnet")) {

      full.spec = full.spec && all(!is.na(nconfig)) && (length(nconfig) == 1)

      if (full.spec)
        full.spec = full.spec && (nconfig == nfitted)

    }#THEN

    # do not trigger a warning if no residuals or fitted values are specified.
    if (!full.spec && any(!is.na(nresid) | !is.na(nfitted)))
      warning("different nodes have different number of residuals or fitted values, disregarding.")

  }#THEN

  # cross-check distributions for consistency and populate the bn.fit object.
  for (node in nodes) {

    if (is(fitted[[node]], c("bn.fit.dnode", "bn.fit.onode"))) {

      # cross-check the levels of each node across all CPTs.
      cpt.levels = lapply(dist, function(x) dimnames(x)[[1]])

      for (cpd in names(dist)[discrete]) {

        # sanity check the new object by comparing it to the old one.
        dist[[cpd]] = check.dnode.vs.spec(dist[[cpd]], old = fitted[[cpd]]$parents,
                        node = cpd, cpt.levels = cpt.levels)
        # all parents of discrete nodes must be discrete nodes themselves.
        if (any(!discrete[fitted[[cpd]]$parents]))
          stop("node ", node, " is discrete but has continuous parents.")
        # store the new CPT in the bn.fit object.
        fitted[[cpd]]$prob = normalize.cpt(dist[[cpd]])

      }#FOR

    }#THEN
    else {

      if (is(fitted[[node]], "bn.fit.cgnode")) {

        # identify discrete and continuous parents and configurations.
        parents = fitted[[node]]$parents
        configs = as.character(seq(from = 0, to = ncol(dist[[node]]$coef) - 1))
        # store the new coefficients and standard deviations.
        fitted[[node]]$coefficients = noattr(dist[[node]]$coef)
        fitted[[node]]$sd = structure(noattr(dist[[node]]$sd), names = configs)
        # identify discrete and continuous parents.
        dparents = as.integer(which(discrete[parents]))
        gparents = as.integer(which(!discrete[parents]))
        fitted[[node]]$dparents = dparents
        fitted[[node]]$gparents = gparents
        # include the levels of the discrete parents
        dlevels = sapply(parents[dparents],
          function(x) dimnames(dist[[x]])[[1]], simplify = FALSE)
        fitted[[node]]$dlevels = dlevels
        # reset columns names for the coefficients and names for sd.
        colnames(fitted[[node]]$coefficients) =
          names(fitted[[node]]$sd) = configs
        # save the configurations of the discrete parents.
        if (full.spec) {

          # check the number of the discrete parents' configurations is right.
          if (prod(sapply(dlevels, length)) != nlevels(dist[[node]]$configs))
            stop("number of discrete parents configurations for node ", node,
              " is ", nlevels(dist[[node]]$configs), " but should be ",
              prod(sapply(dlevels, length)), ".")

          if (any(levels(dist[[node]]$configs) != configs)) {

            # wrong levels, or wrong order: warn and rename.
            warning("remapping levels of the discrete parents configurations ",
              "for node ", node, ".")

            levels(dist[[node]]$configs) = configs

          }#THEN

          fitted[[node]]$configs = noattr(dist[[node]]$configs)

        }#THEN
        else
          fitted[[node]]$configs = factor(NA, levels = seq(from = 0,
            to = prod(sapply(fitted[[node]]$dlevels, length)) - 1L))

      }#THEN
      else if (is(fitted[[node]], "bn.fit.gnode")) {

        # all parents of Gaussian nodes must be continuous nodes.
        if (any(discrete[fitted[[node]]$parents]))
         stop("node ", node,
           " is Gaussian (not conditional Gaussian) but has discrete parents.")

        # store the new coefficients and standard deviations.
        fitted[[node]]$coefficients = noattr(dist[[node]]$coef, ok = "names")
        fitted[[node]]$sd = noattr(dist[[node]]$sd)

      }#THEN

      if (full.spec) {

        fitted[[node]]$residuals = noattr(dist[[node]]$resid,
                                     ok = character(0))
        fitted[[node]]$fitted.values = noattr(dist[[node]]$fitted,
                                         ok = character(0))

      }#THEN
      else {

        fitted[[node]]$residuals = as.numeric(NA)
        fitted[[node]]$fitted.values = as.numeric(NA)

      }#ELSE

    }#ELSE

  }#FOR

  return(structure(fitted, class = c("bn.fit", secondary.class)))

}#CUSTOM.FIT

# return the corect class based on the node classes.
guess.fitted.class = function(fitted) {

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

