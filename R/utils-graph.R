
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

# convert a set of neighbourhoods (markov blankets, neighbours, ...) into the
# corresponding arc set.
sets2arcs = function(sets, nodes) {

  empty = sapply(sets, function(x) {(length(x$nbr) == 0) || is.null(x$nbr) || identical(x$nbr, "")})
  result = do.call(rbind, lapply(nodes[!empty],
               function(x) { cbind(from = x, to = sets[[x]][['nbr']]) }))

  # return an empty matrix if all the neighbourhoods are empty.
  if (is.null(result))
    matrix(character(0), ncol = 2, dimnames = list(c(), c("from", "to")))
  else
    result

}#SETS2ARCS

# get the root/leaf nodes of a graph.
root.leaf.nodes = function(x, leaf = FALSE) {

  .Call(call_root_nodes,
        bn = x,
        check = leaf)

}#ROOT.LEAF.NODES

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

# apply random arc operators to the graph.
perturb.backend = function(network, iter, ops = c("set", "drop", "reverse"),
    nodes, amat, maxp = Inf, whitelist = NULL, blacklist = NULL,
    debug = FALSE) {

  # use a safety copy of the graph.
  new = network
  # remember the nodes whose score has to be recomputed.
  updates = character(0)

  if (debug) {

    if (!is.null(whitelist))
      cat("  > taking the whitelist into account.\n")
    if (!is.null(blacklist))
      cat("  > taking the blacklist into account.\n")

  }#THEN

  # use a 'for' loop instead of a 'repeat' to avoid infinite loops due to the
  # lack of legal operations.
  for (i in seq_len(3 * iter)) {

    # count the parents of each node.
    nparents = colSums(amat)

    to.be.added = arcs.to.be.added(amat, nodes, blacklist = blacklist,
                    nparents = nparents, maxp = maxp)
    to.be.dropped = arcs.to.be.dropped(new$arcs, whitelist = whitelist)
    to.be.reversed = arcs.to.be.reversed(new$arcs, blacklist = blacklist,
                       nparents = nparents, maxp = maxp)

    # no more operation to do.
    if (iter == 0) break

    # choose which arc operator to use.
    op = sample(ops, 1, replace = TRUE)
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

    # update the adjacency matrix, so that has.path() works on the next
    # iteration.
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
structural.hamming.distance = function(learned, true, wlbl = FALSE,
    cpdag = TRUE, debug = FALSE) {

  if (cpdag) {

    learned = cpdag.backend(learned, wlbl = wlbl)
    true = cpdag.backend(true, wlbl = wlbl)

  }#THEN

  .Call(call_shd,
        learned = learned,
        golden = true,
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
    coll = arcs.unique(coll, nodes = nodes)

  }#THEN

  return(coll)

}#COLLIDERS.BACKEND

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
    tp = how.many(which.tp.dir) + how.many(which.tp.und)/2
    fn = how.many(which.fn.dir) + how.many(which.fn.und)/2
    fp = how.many(which.fp.dir) + how.many(which.fp.und)/2

  }#ELSE

  return(list(tp = tp, fp = fp, fn = fn))

}#COMPARE.BACKEND

# check whether connected components are chordal.
connected.components = function(x, debug = FALSE) {

  # check whether igraph is loaded.
  check.and.load.package("igraph")

  # get the undirected arcs.
  x$arcs = undirected.arcs(x)
  x$nodes = cache.structure(names(x$nodes), arcs = x$arcs)

  # find the connected components...
  connected.components =
    .Call(call_connected_components,
          x = x,
          debug = debug)

  # ... exclude those with just one node, which are automatically chordal...
  relevant.connected.components =
    connected.components[sapply(connected.components, length) > 1]
  # ... and check the rest.
  chordal = sapply(relevant.connected.components, function(nodes) {

    individual.component = subgraph.backend(x, nodes)
    igraph::is_chordal(as.igraph(individual.component))$chordal

  })

  return(list(components = connected.components, chordal = chordal))

}#CONNECTED.COMPONENTS

