
# a modified depth-first search, which is able to cope with cycles
# and mixed/undirected graphs.
has.path = function(from, to, nodes, amat, exclude.direct = FALSE,
    underlying.graph = FALSE, debug = FALSE) {

  .Call(call_has_pdag_path,
        from = which(nodes == from),
        to = which(nodes == to),
        amat = amat,
        nrows = nrow(amat),
        nodes = nodes,
        underlying = underlying.graph,
        exclude.direct = exclude.direct,
        debug = debug)

}#HAS.PATH

# convert a set of neighbourhoods into the corresponding arc set.
mb2arcs = function(mb, nodes) {

  empty.mb = sapply(mb, function(x) {(length(x$nbr) == 0) || is.null(x$nbr) || identical(x$nbr, "")})
  result = do.call(rbind, lapply(nodes[!empty.mb],
               function(x) { cbind(from = x, to = mb[[x]][['nbr']]) }))

  # return an empty matrix if all markov blankets are empty.
  if (is.null(result))
    matrix(character(0), ncol = 2, dimnames = list(c(), c("from", "to")))
  else
    result

}#MB2ARCS

# get the root/leaf nodes of a graph.
root.leaf.nodes = function(x, leaf = FALSE) {

  .Call(call_root_nodes,
        bn = x,
        check = leaf)

}#ROOT.NODES.BACKEND

# backend of mb() for bn.fit objects.
mb.fitted = function(x, node) {

  .Call(call_fitted_mb,
        bn = x,
        target = node)

}#MB.FITTED

# return the skeleton of a graph.
dag2ug.backend = function(x, moral = FALSE, debug = FALSE) {

  nodes = names(x$nodes)

  arcs = .Call(call_dag2ug,
               bn = x,
               moral = moral,
               debug = debug)

  # update the arcs of the network.
  x$arcs = arcs
  # update the network structure.
  x$nodes = cache.structure(nodes = nodes, arcs = arcs)
  # make it really clear the graph does not contain any information about arc
  # directions.
  x$learning$undirected = TRUE

  return(x)

}#DAG2UG.BACKEND

# return a complete orientation of a graph.
pdag2dag.backend = function(arcs, ordering) {

  arcs = .Call(call_pdag2dag,
               arcs = arcs,
               nodes = ordering)

  # check that the new graph is still acyclic.
  if (!is.acyclic(arcs = arcs, nodes = ordering))
    stop("this complete orientation of the graph is not acyclic.")

  return(arcs)

}#PDAG2DAG.BACKEND

# mutilated network graph used in likelihood weighting.
mutilated.backend.bn = function(x, evidence) {

  # this is basically a NOP.
  if (identical(evidence, TRUE))
    return(x)

  # only the node names are used here.
  fixed = names(evidence)
  nodes = names(x$nodes)
  # remove all parents of nodes in the evidence.
  x$arcs = x$arcs[x$arcs[, "to"] %!in% fixed, ]
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

# apply random arc operators to the graph.
perturb.backend = function(network, iter, nodes, amat, whitelist,
    maxp = Inf, blacklist, debug = FALSE) {

  # use a safety copy of the graph.
  new = network
  # remember the nodes whose score has to be recomputed.
  updates = character(0)

  # use a 'for' loop instead of a 'repeat' to avoid the threat of
  # infinite loops (due to the lack of legal operations).
  for (i in seq_len(3 * iter)) {

    # count the parents of each node.
    nparents = colSums(amat)

    to.be.added = arcs.to.be.added(amat, nodes, blacklist = blacklist,
                    nparents = nparents, maxp = maxp)
    to.be.dropped = arcs.to.be.dropped(new$arcs, whitelist)
    to.be.reversed = arcs.to.be.reversed(new$arcs, blacklist, nparents, maxp)

    # no more operation to do.
    if (iter == 0) break

    # choose which arc operator to use.
    op = sample(c("set", "drop", "reverse"), 1, replace = TRUE)
    # choose the arc we apply the operator to.
    if ((op == "set") && (nrow(to.be.added) > 0)) {

      a = sample(seq_len(nrow(to.be.added)), 1, replace = TRUE)
      arc = to.be.added[a, ]

      # if the arc creates cycles, choose another one.
      if (!has.path(arc[2], arc[1], nodes, amat)) {

        new$arcs = set.arc.direction(from = arc[1], to = arc[2],
                     arcs = new$arcs, debug = debug)

        updates = c(updates, arc[2])

      }#THEN

    }#THEN
    else if ((op == "drop") && (nrow(to.be.dropped) > 0)) {

      a = sample(seq_len(nrow(to.be.dropped)), 1, replace = TRUE)
      arc = to.be.dropped[a, ]

      new$arcs = drop.arc.backend(arcs = new$arcs, dropped = arc,
                   debug = debug)

      updates = c(updates, arc[2])

    }#THEN
    else if ((op == "reverse") && (nrow(to.be.reversed) > 0)) {

      a = sample(seq_len(nrow(to.be.reversed)), 1, replace = TRUE)
      arc = to.be.reversed[a, ]

      if (!has.path(arc[1], arc[2], nodes, amat, exclude.direct = TRUE)) {

        new$arcs = reverse.arc.backend(from = arc[1], to = arc[2],
                     arcs = new$arcs, debug = debug)

        updates = c(updates, arc)

      }#THEN

    }#ELSE

    # update the adjacency matrix, so that has.path() works on the next iteration.
    amat = arcs2amat(arcs = new$arcs, nodes = nodes)

    # decrease the iteration counter.
    if (!identical(new$arcs, network$arcs))
      iter = iter - 1

  }#FOR

  # save the names of the nodes whose score is to be updated.
  updates = unique(updates)
  new$updates = array(rep(0, length(updates)), dimnames = list(updates))

  return(new)

}#PERTURB.BACKEND

# structural hamming distance backend.
structural.hamming.distance = function(learned, true, wlbl = FALSE, debug = FALSE) {

  .Call(call_shd,
        learned = cpdag.backend(learned, wlbl = wlbl),
        golden = cpdag.backend(true, wlbl = wlbl),
        debug = debug)

}#STRUCTURAL.HAMMING.DISTANCE

# hamming distance backend.
hamming.distance = function(learned, true, debug = FALSE) {

  .Call(call_shd,
        learned = dag2ug.backend(learned),
        golden = dag2ug.backend(true),
        debug = debug)

}#HAMMING.DISTANCE

# backend for extracting colliders from a network.
colliders.backend = function(x, return.arcs = FALSE, including.shielded = TRUE,
    including.unshielded = TRUE, debug = FALSE) {

  nodes = names(x$nodes)

  coll = .Call(call_colliders,
               arcs = x$arcs,
               nodes = nodes,
               return.arcs = return.arcs,
               shielded = including.shielded,
               unshielded = including.unshielded,
               debug = debug)

  if (return.arcs) {

    coll = arcs.rbind(coll[, c("X", "Z")], coll[, c("Y", "Z")])
    coll = unique.arcs(coll, nodes = nodes)

  }#THEN

  return(coll)

}#COLLIDERS.BACKEND

# test equality of two graphs (nodes and arcs).
equal.backend.bn = function(target, current) {

  .Call(call_all_equal_bn,
        target = target,
        current = current)

}#EQUAL.BACKEND.BN

# test equality of two fitted networks (structure and parameters).
equal.backend.fit = function(target, current, tolerance) {

  # check whether the networks have the same structure.
  same.structure = all.equal(bn.net(target), bn.net(current))

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

      # sanity check the target distribuution by comparing it to the old one.
      tprob = check.dnode.vs.dnode(tprob, cnode)

      # checking that the conditional probability tables are identical.
      if (!isTRUE(all.equal(tprob, cprob, tolerance = tolerance)))
        return(paste("Different probabilities for node", node))

    }#THEN
    else if (target.type %in% c("bn.fit.gnode", "bn.fit.cgnode")) {

      tparams = list(coef = tnode$coefficients, sd = tnode$sd,
                  dlevels = tnode$dlevels)

      if (target.type == "bn.fit.gnode")
        tparams = check.gnode.vs.gnode(tparams, cnode)
      if (target.type == "bn.fit.cgnode")
        tparams = check.cgnode.vs.cgnode(tparams, cnode)

      # checking that the regression coefficients are identical.
      if (!isTRUE(all.equal(tparams$coef, cnode$coefficients, tolerance = tolerance)))
        return(paste("Different regression coefficients for node", node))
      # checking that the standard deviations are identical.
      if (!isTRUE(all.equal(tparams$sd, cnode$sd, tolerance = tolerance)))
        return(paste("Different standard errors for node", node))

      # do not check fitted values, residuals and configurations, or networks
      # fitted from different data sets will never be considered equal.

    }#THEN

  }#FOR

  return(TRUE)

}#EQUAL.BACKEND.FIT

# blacklists based on tiers and orderings.
tiers.backend = function(nodes, debug = FALSE) {

  .Call(call_tiers,
        nodes = nodes,
        debug = debug)

}#TIERS.BACKEND

# generate a subgraph spanning a subset of nodes.
subgraph.backend = function(x, nodes, preserve.learning = FALSE) {

  spanned.arcs = function(x) all(!is.na(match(x, nodes)))

  # create the subgraph spanning the subset of nodes.
  res = empty.graph.backend(nodes)
  # identify which arcs are part of the subgraph.
  spanning = apply(x$arcs, 1, spanned.arcs)
  # update the arcs of the subgraph.
  res$arcs = x$arcs[spanning, , drop = FALSE]
  # update the network structure.
  res$nodes = cache.structure(nodes, arcs = res$arcs)

  if (preserve.learning) {

    # copy all the metadata.
    res$learning = x$learning

    # remove arcs incident on nodes that are no longer there from whitelists,
    # blacklists, illegal arcs in a conditional Gaussian BN.
    for (el in c("whitelist", "blacklist", "illegal")) {

      if (is.null(x$learning[[el]]))
        next

      spanning = apply(x$learning[[el]], 1, spanned.arcs)
      res$learning[[el]] = x$learning[[el]][spanning, , drop = FALSE]

    }#FOR

    # remove arcs incident on nodes that are no longer there from Castelo &
    # Siebes prior.
    if (!is.null(x$learning$args$prior) && (x$learning$args$prior == "cs")) {

      spanning = apply(x$learning$args$beta[, c("from", "to")], 1, spanned.arcs)
      res$learning$args$beta = res$learning$args$beta[spanning, , drop = FALSE]

    }#THEN

  }#THEN

  return(res)

}#SUBGRAPH.BACKEND

# test d-separation.
dseparation = function(bn, x, y, z) {

  # this function implements the algorithm from Koller & Friedman's
  # "Probabilistic Graphical Models", Sec. 4.5, pages 136-137.

  # if either x or y are in z, conditional independence is satisfied.
  if ((x %in% z) || (y %in% z))
    return(TRUE)

  # construct the upper closure of the query nodes.
  upper.closure = topological.ordering(bn, start = c(x, y, z), reverse = TRUE)
  ucgraph = subgraph.backend(bn, upper.closure)
  # get its moral graph.
  mgraph = dag2ug.backend(ucgraph, moral = TRUE)
  # remove z and incident arcs to block paths.
  upper.closure = upper.closure[upper.closure %!in% z]
  mgraph = subgraph.backend(mgraph, upper.closure)
  # look for a path between x and y.
  connected = has.path(x, y, nodes = upper.closure,
                amat = arcs2amat(mgraph$arcs, upper.closure))

  return(!connected)

}#DSEPARATION

# compare the arcs in two graphs.
compare.backend = function(target.arcs, current.arcs, nodes, arcs = FALSE) {

  # separate directed and undirected arcs in the target network.
  which.dir = which.directed(target.arcs, nodes)
  target.dir = target.arcs[which.dir, , drop = FALSE]
  target.und = target.arcs[!which.dir, , drop = FALSE]
  # separate directed and undirected arcs in the current network.
  which.dir = which.directed(current.arcs, nodes)
  current.dir = current.arcs[which.dir, , drop = FALSE]
  current.und = current.arcs[!which.dir, , drop = FALSE]

  # treat directed and undirected arcs separately; directed arcs are
  # compared in both presence and direction.
  which.tp.dir = which.listed(target.dir, current.dir)
  which.tp.und = which.listed(target.und, current.und)
  which.fn.dir = !which.listed(target.dir, current.dir)
  which.fn.und = !which.listed(target.und, current.und)
  which.fp.dir = !which.listed(current.dir, target.dir)
  which.fp.und = !which.listed(current.und, target.und)

  if (arcs) {

    # return the arcs corresponding to each category.
    tp = arcs.rbind(target.dir[which.tp.dir, , drop = FALSE],
                    target.und[which.tp.und, , drop = FALSE])
    fn = arcs.rbind(target.dir[which.fn.dir, , drop = FALSE],
                    target.und[which.fn.und, , drop = FALSE])
    fp = arcs.rbind(current.dir[which.fp.dir, , drop = FALSE],
                    current.und[which.fp.und, , drop = FALSE])

  }#THEN
  else {

    # return the counts for each category.
    tp = length(which(which.tp.dir)) + length(which(which.tp.und))/2
    fn = length(which(which.fn.dir)) + length(which(which.fn.und))/2
    fp = length(which(which.fp.dir)) + length(which(which.fp.und))/2

  }#ELSE

  return(list(tp = tp, fp = fp, fn = fn))

}#COMPARE.BACKEND

