
# check one discrete local distribution against another.
check.dnode.vs.dnode = function(new, old) {

  # check that the two sets of parents match, and reorder to match that in the
  # old network.
  if (!setequal(new$parents, old$parents))
    stop("the node", new$node, "has different parents in the two networks.")

  new$parents = old$parents

  # do not check children at all, just replace them with those in the old
  # network.
  new$children = old$children

  # check and rearrange the conditional probability table.
  new$prob = check.rvalue.vs.dnode(new$prob, old)

  return(new)

}#CHECK.DNODE.VS.DNODE

# check one Gaussian local distribution against another.
check.gnode.vs.gnode = function(new, old, keep.fitted = FALSE) {

  # check that the two sets of parents match, and reorder to match that in the
  # old network.
  if (!setequal(new$parents, old$parents))
    stop("the node", new$node, "has different parents in the two networks.")

  new$parents = old$parents

  # do not check children at all, just replace them with those in the old
  # network.
  new$children = old$children

  if (keep.fitted) {

    rvalue = list(coef = new$coefficients, sd = new$sd,
                  fitted = new$fitted.values, resid = new$residuals)

  }#THEN
  else {

    rvalue = list(coef = new$coefficients, sd = new$sd,
                  fitted = NULL, resid = NULL)

  }#ELSE

  rvalue = check.rvalue.vs.gnode(rvalue, old)

  new$coefficients = rvalue$coef
  new$sd = rvalue$sd
  new$fitted.values = rvalue$fitted
  new$residuals = rvalue$residuals

  return(new)

}#CHECK.GNODE.VS.GNODE

check.cgnode.vs.cgnode = function(new, old, keep.fitted = FALSE) {

  # check that the two sets of parents match, and reorder to match that in the
  # old network.
  if (!setequal(new$parents, old$parents))
    stop("the node", new$node, "has different parents in the two networks.")

  new$parents = old$parents

  # do not check children at all, just replace them with those in the old
  # network.
  new$children = old$children

  if (keep.fitted) {

    rvalue = list(coef = new$coefficients, sd = new$sd,
                  fitted = new$fitted.values, resid = new$residuals,
                  configs = new$configs)

  }#THEN
  else {

    rvalue = list(coef = new$coefficients, sd = new$sd,
                  fitted = NULL, resid = NULL, configs = NULL)

  }#ELSE

  rvalue = check.rvalue.vs.cgnode(rvalue, old)

  new$coefficients = rvalue$coef
  new$sd = rvalue$sd
  new$fitted.values = rvalue$fitted
  new$residuals = rvalue$residuals
  new$configs = rvalue$configs

  return(new)

}#CHECK.CGNODE.VS.CGNODE

# does the network has any NA parameters?
is.bn.fit.ill.defined = function(fitted) {

  for (node in names(fitted)) {

    fit = fitted[[node]]

    if (is(fit, c("bn.fit.dnode", "bn.fit.onode"))) {

      if (anyNA(fit$prob))
        return(TRUE)

    }#THEN
    else if (is(fit, c("bn.fit.gnode", "bn.fit.cgnode"))) {

      if (anyNA(fit$coefficients))
        return(TRUE)
      if (anyNA(fit$sd))
        return(TRUE)

    }#THEN

  }#FOR

  return(FALSE)

}#IS.BN.FIT.ILL.DEFINED

# check one bn.fit object against another.
check.fitted.vs.fitted = function(new, old, local = TRUE) {

  # they should have the same class.
  class.new = grep("net$", class(new), value = TRUE)
  class.old = grep("net$", class(old), value = TRUE)
  if (class.new != class.old)
    stop("the first bn.fit object has class '", class.new,
         "', the second has class '", class.old, "'.")

  # compare two networks at the level of the local distributions.
  if (local) {

    # they should have the same structure.
    if (!isTRUE(all.equal(bn.net(new), bn.net(old))))
      stop("the two bn.fit objects have different underlying graphs.")

    all.classes.new = class(new)
    class(new) = "list"

    # the nodes should have the same parameterisation.
    for (node in names(new)) {

      node.new = new[[node]]
      node.old = old[[node]]

      if (class(node.new) != class(node.new))
        stop("node", node, "is a '", class(node.new), "' in P and '",
             class(node.new), "' in Q.")

      if (is(node.new, c("bn.fit.dnode", "bn.fit.onode")))
        node.new = check.dnode.vs.dnode(node.new, node.old)
      else if (is(node.new, "bn.fit.gnode"))
        node.new = check.gnode.vs.gnode(node.new, node.old)
      else if (is(node.new, "bn.fit.cgnode"))
        node.new = check.cgnode.vs.cgnode(node.new, node.old)

      new[[node]] = node.new

    }#FOR

    # store the nodes in the same order in both networks.
    new = new[names(old)]

    class(new) = all.classes.new

  }#THEN
  else {

    all.classes.new = class(new)
    class(new) = "list"

    # compare two networks at the level of the global distribution: they should
    # have the same nodes...
    if (!setequal(names(new), names(old)))
      stop("the two bn.fit objects have different nodes.")
    # ... in the same order...
    new = new[names(old)]
    # ... each node...
    for (node in names(new)) {

      node.new = new[[node]]
      node.old = old[[node]]

      # ... should be either a continuous or a discrete variables in both
      # networks...
      discrete.new = is(node.new, c("bn.fit.dnode", "bn.fit.onode"))
      discrete.old = is(node.old, c("bn.fit.dnode", "bn.fit.onode"))

      if (discrete.new && !discrete.old)
        stop("node ", node, " is discrete in fhe first bn.fit object but not in the second.")
      if (!discrete.new && discrete.old)
        stop("node ", node, " is discrete in fhe second bn.fit object but not in the first.")

      # ... and if they are both discrete they should have the same levels.
      if (discrete.new && discrete.old) {

        levels.new = dimnames(node.new$prob)[[1]]
        levels.old = dimnames(node.old$prob)[[1]]

        if (!setequal(levels.new, levels.old))
          stop("node ", node, "has different levels in the two bn.fit objects.")

      }#THEN

    }#FOR

    class(new) = all.classes.new

  }#THEN

  return(new)

}#CHECK.FITTED.VS.FITTED

