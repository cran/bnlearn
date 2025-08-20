# mutilated network graph used in likelihood weighting.
mutilated.backend.bn = function(x, evidence) {

  # this is basically a NOP.
  if (identical(evidence, TRUE))
    return(x)

  # only the node names are used here.
  fixed = names(evidence)
  nodes = names(x$nodes)
  # remove all parents of nodes in the evidence.
  x$arcs = x$arcs[x$arcs[, "to"] %!in% fixed, , drop = FALSE]
  # update the cached information for the fixed nodes.
  amat = arcs2amat(x$arcs, nodes)
  for (node in fixed)
    x$nodes[[node]] = cache.partial.structure(nodes, target = node,
                        amat = amat, debug = FALSE)

  return(x)

}#MUTILATED.BACKEND.BN

# mutilated fitted network used in likelihood sampling.
mutilated.backend.fitted = function(x, evidence) {

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
      if (!is.null(cur$fitted.values))
        cur$residuals = rep(fix, length(cur$fitted.values))
      if (!is.null(cur$residuals))
        cur$residuals = rep(0, length(cur$residuals))

    }#THEN
    else if (is(cur, c("bn.fit.dnode", "bn.fit.onode"))) {

      # reset the conditional distribution.
      levels = dimnames(cur$prob)[[1]]
      values = (levels %in% fix) + 0

      cur$prob = as.table(structure(values / sum(values), names = levels))

    }#THEN

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

}#MUTILATED.BACKEND.FITTED

# twin network used for counterfactuals.
twin.backend.bn = function(x) {

  # function to map factual nodes to counterfactual and exogenous nodes.
  ctf = function(nodes) paste0(nodes, ".", recycle0 = TRUE)
  exo = function(nodes) paste0("u", nodes, recycle0 = TRUE)
  # function to map counterfactual to factual nodes.
  fac = function(nodes) sub("\\.$", "", nodes)

  # extract the node and arc set of the original network...
  orig.nodes = .nodes(x)
  orig.arcs = x$arcs
  # ... create the counterfactual and exogenous nodes...
  counterfactual.nodes = ctf(orig.nodes)
  exogenous.nodes = exo(orig.nodes)
  all.nodes = c(orig.nodes, counterfactual.nodes, exogenous.nodes)
  # ... verify that node labels are still unique...
  dupes = duplicated(all.nodes)
  if (any(dupes))
    stop("duplicated node labels in twin network: ",
         paste(paste0("'", unique(all.nodes[dupes]), "'"), collapse = ", "), ".")

  # ... replicate the original arcs for the counterfactual nodes...
  counterfactual.arcs =
    matrix(ctf(x$arcs), ncol = 2, dimnames = list(NULL, c("from", "to")))
  # connect the exogenous nodes to the corresponding {counter,}factual nodes.
  exogenous.arcs =
    matrix(c(rep(exogenous.nodes, 2), c(orig.nodes, counterfactual.nodes)),
           ncol = 2, dimnames = list(NULL, c("from", "to")))

  all.arcs = rbind(orig.arcs, counterfactual.arcs, exogenous.arcs)

  # create the bn object.
  twin.network = empty.graph.backend(all.nodes)
  twin.network$arcs = all.arcs
  twin.network$nodes = cache.structure(all.nodes, all.arcs)

  # set the class, and remove all other optional classes (to avoid misapplying
  # methods for classifiers, for instance).
  class(twin.network) = c("bn", "bn.twin")
  # set additional attributes for printing, interventions and node merging.
  twin.network$learning$algo = "twin"
  twin.network$learning$label.mapping = c(ctf = ctf, exo = exo, fac = fac)
  twin.network$learning$causal.roles =
    list(factual = orig.nodes, exogenous = exogenous.nodes,
      counterfactual = counterfactual.nodes)

  return(twin.network)

}#TWIN.BACKEND.BN

# set up a counterfactual network.
counterfactual.backend.bn = function(x, evidence, merging = TRUE) {

  # create the full twin network.
  full.twin = twin.backend.bn(x)
  # perform the intervention to the counterfactual node correspoding to the node
  # in the evidence.
  map = full.twin$learning$label.mapping
  names(evidence) = map$ctf(names(evidence))
  ctf.twin = mutilated.backend.bn(full.twin, evidence = evidence)

  if (!merging)
    return(ctf.twin)

  # merge nodes that have the same parents, equivalently, drop counterfactual
  # node that have the same distribution as their factual counterparts and
  # rewire ourgoing arcs.
  ctf.nodes = ctf.twin$learning$causal.roles$counterfactual
  to.drop = character(0)

  for (node in ctf.nodes) {

    # counterfactual nodes that are target of interventions are now different
    # from the corresponding factual node.
    if (node %in% names(evidence))
      next

    # therefore, the children of such nodes cannot be merged either.
    if (any(ctf.twin$nodes[[node]]$parents %in% names(evidence)))
      next

    # for other nodes, remove both the exogenous and the counterfactual node.
    to.drop = c(to.drop, node)

  }#FOR

  for (node in ctf.nodes) {

    merge.nodes = intersect(ctf.twin$nodes[[node]]$parents, to.drop)
    for (merge in merge.nodes)
      ctf.twin = set.arc(ctf.twin, from = map$fac(merge), to = node)

  }#FOR

  reduced.twin = subgraph.backend(ctf.twin, setdiff(nodes(ctf.twin), to.drop))

  return(reduced.twin)

}#COUNTERFACTUAL.BACKEND.BN
