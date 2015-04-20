
# a modified depth-first search, which is able to cope with cycles
# and mixed/undirected graphs.
has.path = function(from, to, nodes, amat, exclude.direct = FALSE,
    underlying.graph = FALSE, debug = FALSE) {

  .Call("has_pdag_path",
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

  .Call("root_nodes",
        bn = x,
        check = as.integer(leaf))

}#ROOT.NODES.BACKEND

# get the parents of a node.
parents.backend = function(arcs, node, undirected = FALSE) {

  if (!undirected)
    arcs[(arcs[, "to"] == node) & which.directed(arcs), "from"]
  else
    arcs[(arcs[, "to"] == node), "from"]

}#PARENTS.BACKEND

# backend of mb() for bn.fit objects.
mb.fitted = function(x, node) {

  .Call("fitted_mb",
        bn = x,
        target = node)

}#MB.FITTED

# return the skeleton of a graph.
dag2ug.backend = function(x, moral = FALSE, debug = FALSE) {

  nodes = names(x$nodes)

  arcs = .Call("dag2ug",
               bn = x,
               moral = moral,
               debug = debug)

  # update the arcs of the network.
  x$arcs = arcs
  # update the network structure.
  x$nodes = cache.structure(nodes = nodes, arcs = arcs)

  return(x)

}#DAG2UG.BACKEND

# return a complete orientation of a graph.
pdag2dag.backend = function(arcs, ordering) {

  arcs = .Call("pdag2dag",
               arcs = arcs,
               nodes = ordering)

  # check that the new graph is still acyclic.
  if (!is.acyclic(arcs = arcs, nodes = ordering))
    stop("this complete orientation of the graph is not acyclic.")

  return(arcs)

}#PDAG2DAG.BACKEND

# reconstruct the equivalence class of a network.
cpdag.backend = function(x, moral = TRUE, debug = FALSE) {

  nodes = names(x$nodes)

  amat = .Call("cpdag",
               arcs = x$arcs,
               nodes = nodes,
               moral = moral,
               fix = FALSE,
               debug = debug)

  # update the arcs of the network.
  x$arcs = amat2arcs(amat, nodes)

  # update the network structure.
  x$nodes = cache.structure(nodes, amat = amat)

  return(x)

}#CPDAG.BACKEND

# reconstruct the arc set of the equivalence class of a network.
cpdag.arc.backend = function(nodes, arcs, moral = TRUE, fix.directed = FALSE,
    debug = FALSE) {

  amat = .Call("cpdag",
               arcs = arcs,
               nodes = nodes,
               moral = moral,
               fix = fix.directed,
               debug = debug)

  return(amat2arcs(amat, nodes))

}#CPDAG.ARCS.BACKEND

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
    else if(is(cur, c("bn.fit.dnode", "bn.fit.onode"))) {

      # reset the conditional distribution.
      cur$prob = as.table(structure((dimnames(cur$prob)[[1]] == fix) + 0,
                   names = dimnames(cur$prob)[[1]]))

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
structural.hamming.distance = function(learned, true, debug = FALSE) {

  .Call("shd",
        learned = cpdag.backend(learned),
        golden = cpdag.backend(true),
        debug = debug)

}#STRUCTURAL.HAMMING.DISTANCE

# hamming distance backend.
hamming.distance = function(learned, true, debug = FALSE) {

  .Call("shd",
        learned = dag2ug.backend(learned),
        golden = dag2ug.backend(true),
        debug = debug)

}#HAMMING.DISTANCE

# backend for extracting v-structures from a network.
vstructures = function(x, arcs, moral = TRUE, debug = FALSE) {

  .Call("vstructures",
        arcs = x$arcs,
        nodes = names(x$nodes),
        return.arcs = arcs,
        moral = moral,
        debug = debug)

}#VSTRUCTURES

# test equality of two graphs (nodes and arcs).
equal.backend = function(target, current) {

  .Call("all_equal",
        target = target,
        current = current)

}#EQUAL.BACKEND

# blacklists based on tiers and orderings.
tiers.backend = function(nodes, debug = FALSE) {

  .Call("tiers",
        nodes = nodes,
        debug = debug)

}#TIERS.BACKEND

# backend to get a DAG out of a CPDAG (still in the same equivalence class).
cpdag.extension = function(x, debug = FALSE) {

  nodes = names(x$nodes)

  # update the arcs of the network.
  x$arcs = cpdag.arc.extension(arcs = x$arcs, nodes = nodes, debug = debug)
  # update the network structure.
  x$nodes = cache.structure(nodes, arcs = x$arcs, debug = debug)

  return(x)

}#CPDAG.EXTENSION

# backend to get a set of directed arcs out of a CPDAG.
cpdag.arc.extension = function(arcs, nodes, debug = FALSE) {

  .Call("pdag_extension",
        arcs = arcs,
        nodes = nodes,
        debug = debug)

}#CPDAG.ARC.EXTENSION

# generate a subgraph spanning a subset of nodes.
subgraph.backend = function(x, nodes) {

  # create the subgraph spanning the subset of nodes.
  res = empty.graph.backend(nodes)
  # identify which arcs are part of the subgraph.
  spanning = apply(x$arcs, 1, function(x) all(!is.na(match(x, nodes))))
  # update the arcs of the subgraph.
  res$arcs = x$arcs[spanning, , drop = FALSE]
  # update the network structure.
  res$nodes = cache.structure(nodes, arcs = res$arcs)

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
  upper.closure = schedule(bn, start = c(x, y, z), reverse = TRUE)
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
