# intervention on a structural causal model.
intervention.backend = function(scm, evidence) {

  # this is basically a NOP.
  if (identical(evidence, TRUE))
    return(scm)

  # only the node names are used here.
  fixed = names(evidence)
  nodes = names(scm$nodes)
  # remove all parents of nodes in the evidence.
  scm$arcs = scm$arcs[scm$arcs[, "to"] %!in% fixed, , drop = FALSE]
  # update the metadata.
  for (f in fixed) {

    # remove the link to the exogenous parent.
    exo = scm$nodes[[f]]$exogenous
    scm$nodes[[f]]$exogenous = character(0)
    scm$nodes[[exo]]$factual = character(0)
    # remove the link to the factual parents.
    for (pa in scm$nodes[[f]]$parents) {

      scm$nodes[[pa]]$children = setdiff(scm$nodes[[pa]]$children, f)
      scm$nodes[[f]]$parents = character(0)

    }#FOR

  }#FOR

  return(scm)

}#INTERVENTION.BACKEND

# mutilated fitted network used for causal inference and likelihood sampling.
intervention.backend.fitted = function(x, evidence) {

  # this is basically a NOP.
  if (identical(evidence, TRUE))
    return(x)

  # extract the names of the nodes.
  fixed = names(evidence)
  nodes = names(x)

  for (node in fixed) {

    # cache the node information.
    cur = x[[node]]
    fix = evidence[[node]]

    if (is(cur, "bn.fit.gnode")) {

      # reset the conditional distribution.
      cur$coefficients = c("(Intercept)" = as.numeric(fix))
      cur$sd = 0
      # reset fitted values and residuals.
      if (!is.null(cur$fitted.values) && !all(is.na(cur$fitted.values)))
        cur$residuals = rep(fix, length(cur$fitted.values))
      if (!is.null(cur$residuals) && !all(is.na(cur$residuals)))
        cur$residuals = rep(0, length(cur$residuals))

    }#THEN
    else if (is(cur, c("bn.fit.dnode", "bn.fit.onode"))) {

      # reset the conditional distribution.
      levels = dimnames(cur$prob)[[1]]
      values = (levels %in% fix) + 0

      cur$prob = as.table(structure(values / sum(values), names = levels))

    }#THEN
    else if (is(cur, "bn.fit.cgnode")) {

      # reset a conditional gaussian distribution to that of a gaussian node.
      class(cur) = "bn.fit.gnode"
      cur$dparents = cur$gparents = cur$configs = cur$dlevels = NULL
      cur$coefficients = c("(Intercept)" = as.numeric(fix))
      cur$sd = 0
      # reset fitted values and residuals.
      if (!is.null(cur$fitted.values) && !all(is.na(cur$fitted.values)))
        cur$residuals = rep(fix, length(cur$fitted.values))
      if (!is.null(cur$residuals) && !all(is.na(cur$residuals)))
        cur$residuals = rep(0, length(cur$residuals))

    }#THEN
    else if (is(cur, "bn.fit.zihpnode")) {

      if (fix == 0) {

        # structural zeroes with probability 1.
        cur$inflation = c("(Intercept)" = Inf)
        cur$intensity = c("(Intercept)" = 1)
        cur$dispersion = 1

      }#THEN
      else {

        # no zero inflation, hyper-Poisson fixed to a single value.
        cur$inflation = c("(Intercept)" = -Inf)
        cur$intensity = c("(Intercept)" = log(as.numeric(fix)))
        cur$dispersion = -Inf

      }#THEN

    }#THEN
    else if (is(cur, "bn.fit.zinbnode")) {

      if (fix == 0) {

        # structural zeroes with probability 1.
        cur$inflation = c("(Intercept)" = Inf)
        cur$succprob = c("(Intercept)" = 0)
        cur$failures = 0

      }#THEN
      else {

        stop("unsupported.")

      }#THEN

    }#THEN
    else
      stop("node ", node, " has unknown class ", class(cur), ".")

    # update parents and children.
    parents = cur$parents
    cur$parents = character(0)

    for (p in parents) {

      temp = x[[p]]
      temp$children = temp$children[temp$children != node]
      x[p] = list(temp)

    }#FOR

    x[node] = list(cur)

  }#FOR

  return(x)

}#INTERVENTION.BACKEND.FITTED

# twin network for a structural causal model.
twin.backend = function(scm) {

  # already a twin network, nothing to do.
  if (is(scm, "scm.twin"))
    return(scm)

  # add counterfactual node labels.
  scm$roles$counterfactual = .ctf(scm$roles$factual)
  # verify that node labels are still unique.
  all.nodes = unlist(scm$roles)
  dupes = duplicated(all.nodes)
  if (any(dupes))
    stop("duplicated node labels in twin network: ",
         paste(paste0("'", unique(all.nodes[dupes]), "'"), collapse = ", "), ".")

  # add the metadata for the counterfactual nodes.
  for (f in scm$roles$factual)
    scm$nodes[[.ctf(f)]] = list(
      counterfactual = character(0),
      exogenous = scm$nodes[[f]]$exogenous,
      factual = f,
      parents = .ctf(scm$nodes[[f]]$parents),
      children = .ctf(scm$nodes[[f]]$children)
    )
  # link the factual nodes to the counterfactual nodes.
  for (f in scm$roles$factual)
    scm$nodes[[f]]$counterfactual = .ctf(f)
  # link the exogenous nodes to the counterfactual nodes.
  for (e in scm$roles$exogenous)
    scm$nodes[[e]]$counterfactual = .ctf(scm$nodes[[e]]$factual)

  # duplicate the arcs between factual nodes for the counterfactual nodes.
  selected = (scm$arcs[, "from"] %in% scm$roles$factual) &
             (scm$arcs[, "to"] %in% scm$roles$factual)
  factual.arcs = counterfactual.arcs = scm$arcs[selected, , drop = FALSE]
  counterfactual.arcs[] = .ctf(counterfactual.arcs)
  # duplicate the arcs from the exogenous nodes.
  selected = (scm$arcs[, "from"] %in% scm$roles$exogenous) &
             (scm$arcs[, "to"] %in% scm$roles$factual)
  exo.to.fac.arcs = exo.to.ctf.arcs = scm$arcs[selected, , drop = FALSE]
  exo.to.ctf.arcs[, "to"] = .ctf(exo.to.ctf.arcs[, "to"])

  # collate all these subsets of arcs together.
  scm$arcs = rbind(factual.arcs, counterfactual.arcs,
                   exo.to.fac.arcs, exo.to.ctf.arcs)

  # update the classes, to match the behaviour for bn objects.
  class(scm) = c("scm", "scm.twin")

  return(scm)

}#TWIN.BACKEND

# twin network used for counterfactuals, with parameters.
twin.backend.fitted = function(x) {

  # already a twin network, nothing to do.
  if (is(x, "bn.twin"))
    return(x)

  # extract the node set of the original network and remove its class to avoid
  # method dispatch...
  orig.nodes = .nodes(x)
  x = unclass(x)
  # ... create the counterfactual and exogenous nodes...
  counterfactual.nodes = .ctf(orig.nodes)
  exogenous.nodes = .exo(orig.nodes)
  all.nodes = c(orig.nodes, counterfactual.nodes, exogenous.nodes)
  # ... verify that node labels are still unique.
  dupes = duplicated(all.nodes)
  if (any(dupes))
    stop("duplicated node labels in twin network: ",
         paste(paste0("'", unique(all.nodes[dupes]), "'"), collapse = ", "), ".")

  twin.network = vector(length(all.nodes), mode = "list")
  names(twin.network) = all.nodes

  # for each node in the original network:
  for (node in orig.nodes) {

    # the corresponding factual node has the same structure, without the
    # error term and with the associated exogenous node as an additional
    # parent.
    factual = x[[node]]
    factual$parents = c(factual$parents, .exo(node))
    factual$coefficients[.exo(node)] = 1
    factual$sd = 0
    factual$residuals = factual$fitted.values = NA_real_

    # the corresponding exogenous node is a root node with mean zero and the
    # standard error from the original node, with that node and its
    # counterfactual as children.
    exogenous = list(
      node = .exo(node),
      parents = character(0),
      children = c(node, .ctf(node)),
      coefficients = c("(Intercept)" = 0),
      sd = x[[node]]$sd,
      residuals = NA_real_,
      fitted.values = NA_real_
    )
    class(exogenous) = "bn.fit.gnode"

    # the corresponding counterfactual node is similar to the factual node,
    # with the factual parents/children replaced by the corresponding
    # counterfactual nodes.
    counterfactual = x[[node]]
    counterfactual$node = .ctf(node)
    counterfactual$parents = c(.ctf(counterfactual$parents), .exo(node))
    counterfactual$children = .ctf(counterfactual$children)
    names(counterfactual$coefficients)[-1] =
      .ctf(names(counterfactual$coefficients)[-1])
    counterfactual$coefficients[.exo(node)] = 1
    counterfactual$sd = 0
    counterfactual$residuals = counterfactual$fitted.values = NA_real_

    # always drop the fitted values and the residuals from the original network,
    # treat the twin network as a new, separate model.
    twin.network[[node]] = factual
    twin.network[[.exo(node)]] = exogenous
    twin.network[[.ctf(node)]] = counterfactual

  }#FOR

  # propagate causal roles to facilitate downstream use.
  attr(twin.network, "roles") =
    list(factual = orig.nodes, exogenous = .exo(orig.nodes),
         counterfactual = .ctf(orig.nodes))
  # set the class, and remove all other optional classes (to avoid misapplying
  # methods for classifiers, for instance).
  class(twin.network) = c("bn.fit", "bn.twin", "bn.fit.gnet")

  return(twin.network)

}#TWIN.BACKEND.FITTED

# set up a counterfactual structural causal model.
counterfactual.backend = function(full.twin, evidence, merging = TRUE) {

  # perform the intervention to the counterfactual node corresponding to the node
  # in the evidence.
  ctf.twin = intervention.backend(full.twin, evidence = evidence)
  # perform node merging, if requested.
  if (merging)
    ctf.twin = counterfactual.node.merging(ctf.twin, evidence = evidence)

  class(ctf.twin) = c(class(ctf.twin), "scm.ctf")

  return(ctf.twin)

}#COUNTERFACTUAL.BACKEND

# set up a counterfactual network.
counterfactual.backend.fitted = function(full.twin, evidence, merging = TRUE) {

  # perform the intervention to the counterfactual node corresponding to the node
  # in the evidence.
  ctf.twin = intervention.backend.fitted(full.twin, evidence = evidence)

  # either return the twin network as is...
  if (!merging)
    return(ctf.twin)

  # ... or perform the node merging on the structure.
  twin.dag = bn.net(ctf.twin)
  twin.scm = from.bn.twin.to.scm.twin(twin.dag)
  reduced.scm = counterfactual.node.merging(twin.scm, evidence = evidence)
  reduced.dag = from.scm.to.bn(reduced.scm)
  nodes = .nodes(reduced.dag)

  # copy the parameters and local structure into to the reduced network.
  fitted = structure(vector(length(nodes), mode = "list"), names = nodes)

  for (node in names(fitted)) {

    ldist = ctf.twin[[node]]

    if (!setequal(ldist$parents, reduced.dag$nodes[[node]]$parents)) {

      # some parents are counterfactual nodes that have been merged, replace
      # them with the corresponding factual nodes.
      replace = setdiff(ldist$parents, reduced.dag$nodes[[node]]$parents)
      coef.pos = match(replace, names(ldist$coefficients))
      names(ldist$coefficients)[coef.pos] = .fac.ctf(replace)

    }#ELSE

    # always drop the fitted values and the residuals from the original network,
    # treat the counterfactual network as a new, separate model.
    fitted[[node]] = structure(list(
      node = node,
      parents = reduced.dag$nodes[[node]]$parents,
      children = reduced.dag$nodes[[node]]$children,
      coefficients = ldist$coefficients,
      sd = ldist$sd,
      residuals = NA_real_,
      fitted.values = NA_real_
    ), class = "bn.fit.gnode")

  }#FOR

  return(structure(fitted,
           class = c("bn.fit", "bn.twin", "bn.ctf", "bn.fit.gnet")))

}#COUNTERFACTUAL.BACKEND.FITTED

# optional node merging for twin networks after applying counterfactuals.
counterfactual.node.merging = function(twin, evidence) {

  # merge nodes that have the same parents, equivalently, drop counterfactual
  # node that have the same distribution as their factual counterparts and
  # rewire outgoing arcs.
  ctf.nodes = twin$roles$counterfactual
  to.drop = character(0)

  for (node in ctf.nodes) {

    # counterfactual nodes that are target of interventions are now different
    # from the corresponding factual node.
    if (node %in% names(evidence))
      next

    # therefore, the children of such nodes cannot be merged either.
    if (any(twin$nodes[[node]]$parents %in% names(evidence)))
      next

    # for other nodes, remove the counterfactual node.
    to.drop = c(to.drop, node)

  }#FOR

  # adjust the metadata for each counterfactual node to be dropped by...
  for (node in to.drop) {

    # ... removing it from its factual and exogenous nodes...
    twin$nodes[[ twin$nodes[[node]]$factual ]]$counterfactual = character(0)
    twin$nodes[[ twin$nodes[[node]]$exogenous ]]$counterfactual = character(0)
    # ... reassigning children to the factual node....
    children.to.keep = setdiff(twin$nodes[[node]]$children, to.drop)

    if (length(children.to.keep) > 0) {

      for (chld in children.to.keep) {

        match = which(twin$nodes[[chld]]$parents == node)
        twin$nodes[[chld]]$parents[match] = twin$nodes[[node]]$factual

      }#FOR
      # ... and remapping arcs that point from counterfactual nodes to be
      # dropped to nodes that are not dropped; they should go out from the
      # factual nodes instead.
      remap = (twin$arcs[, "from"] == node) &
              (twin$arcs[, "to"] %in% children.to.keep)
      twin$arcs[remap, "from"] = twin$nodes[[node]]$factual

    }#THEN

    # finally, remove its own metadata.
    twin$nodes[[node]] = NULL

  }#FOR

  # remove all those nodes from the roles.
  twin$roles$counterfactual = setdiff(twin$roles$counterfactual, to.drop)
  # remove all arcs incident on those nodes.
  adrop = (twin$arcs[, "from"] %in% to.drop) | (twin$arcs[, "to"] %in% to.drop)
  twin$arcs = twin$arcs[!adrop, , drop = FALSE]

  return(twin)

}#COUNTERFACTUAL.NODE.MERGING

