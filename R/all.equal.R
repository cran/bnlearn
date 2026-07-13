# test equality of two graphs (nodes and arcs).
equal.backend.bn = function(target, current) {

  .Call(call_all_equal_bn,
        target = target,
        current = current)

}#EQUAL.BACKEND.BN

# test equality of two fitted networks (structure and parameters).
equal.backend.fit = function(target, current, tolerance) {

  # check whether the networks have the same structure.
  same.structure = equal.backend.bn(bn.net(target), bn.net(current))

  if (!isTRUE(same.structure))
    return(same.structure)

  for (node in nodes(target)) {

    # extract the nodes from the structure.
    tnode = target[[node]]
    cnode = current[[node]]
    # check whether the nodes follow the same distribution.
    target.type = class(tnode)
    current.type = class(cnode)

    # check whether the parameters of the local distribution are the same.
    if (target.type != current.type)
      return(paste("Different distributions for node", node))

    if (target.type %in% c("bn.fit.dnode", "bn.fit.onode")) {

      tprob = tnode$prob
      cprob = cnode$prob

      # sanity check the target distribution by comparing it to the old one.
      tprob = check.rvalue.vs.dnode(tprob, cnode)

      # checking that the conditional probability tables are identical.
      if (!isTRUE(all.equal(tprob, cprob, tolerance = tolerance)))
        return(paste("Different probabilities for node", node))

    }#THEN
    else if (target.type %in% c("bn.fit.gnode", "bn.fit.cgnode")) {

      tparams = list(coef = tnode$coefficients, sd = tnode$sd,
                  dlevels = tnode$dlevels)

      if (target.type == "bn.fit.gnode")
        tparams = check.rvalue.vs.gnode(tparams, cnode)
      else if (target.type == "bn.fit.cgnode")
        tparams = check.rvalue.vs.cgnode(tparams, cnode)

      # checking that the regression coefficients are identical.
      if (!isTRUE(all.equal(tparams$coef, cnode$coefficients, tolerance = tolerance)))
        return(paste("Different regression coefficients for node", node))
      # checking that the standard deviations are identical.
      if (!isTRUE(all.equal(tparams$sd, cnode$sd, tolerance = tolerance)))
        return(paste("Different standard errors for node", node))

      # do not check fitted values, residuals and configurations, or networks
      # fitted from different data sets will never be considered equal.

    }#THEN
    else if (target.type == "bn.fit.zihpnode") {

      tparams = tnode[c("inflation", "intensity", "dispersion")]

      # sanity check the target distribution by comparing it to the old one.
      tparams = check.rvalue.vs.zihpnode(tparams, cnode)

      # checking that both sets of regression coefficients are identical.
      if (!isTRUE(all.equal(tparams$inflation, cnode$inflation, tolerance = tolerance)))
        return(paste("Different inflation coefficients for node", node))
      if (!isTRUE(all.equal(tparams$intensity, cnode$intensity, tolerance = tolerance)))
        return(paste("Different intensity coefficients for node", node))
      # checking that the dispersion parameters are identical.
      if (!isTRUE(all.equal(tparams$dispersion, cnode$dispersion, tolerance = tolerance)))
        return(paste("Different dispersion for node", node))

    }#THEN
    else if (target.type == "bn.fit.zinbnode") {

      tparams = tnode[c("inflation", "prsucc", "failures")]

      # sanity check the target distribution by comparing it to the old one.
      tparams = check.rvalue.vs.zinbnode(tparams, cnode)

      # checking that both sets of regression coefficients are identical.
      if (!isTRUE(all.equal(tparams$inflation, cnode$inflation, tolerance = tolerance)))
        return(paste("Different inflation coefficients for node", node))
      if (!isTRUE(all.equal(tparams$prsucc, cnode$prsucc, tolerance = tolerance)))
        return(paste("Different prsucc coefficients for node", node))
      # checking that the intensities are identical.
      if (!isTRUE(all.equal(tparams$failures, cnode$failures, tolerance = tolerance)))
        return(paste("Different numbers of failures for node", node))

    }#THEN

  }#FOR

  return(TRUE)

}#EQUAL.BACKEND.FIT


