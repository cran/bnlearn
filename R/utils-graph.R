
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
        debug = debug,
        PACKAGE = "bnlearn")

}#HAS.PATH

# count the cycles the arc is part of (even in a partially directed graph).
how.many.cycles = function(arc, nodes, amat, debug = FALSE) {

  .Call("how_many_cycles",
        from = which(nodes == arc[2]),
        to = which(nodes == arc[1]),
        amat = amat,
        nrows = nrow(amat),
        nodes = nodes,
        debug = debug,
        PACKAGE = "bnlearn")

}#HOW.MANY.CYCLES

# convert a set of neighbourhoods into the corresponding arc set.
mb2arcs = function(mb, nodes) {

  empty.mb = sapply(mb, function(x) {(length(x$nbr) == 0) || is.null(x$nbr) || identical(x$nbr, "")})
  result = do.call(rbind, lapply(nodes[!empty.mb],
               function(x) { cbind(from = x, to = mb[[x]][['nbr']]) }))

  # return an empty matrix all markov blankets are empty.
  if (is.null(result))
    matrix(character(0), ncol = 2, dimnames = list(c(), c("from", "to")))
  else
    result

}#MB2ARCS

# get the root/leaf nodes of a graph.
root.leaf.nodes = function(x, leaf = FALSE) {

  .Call("root_nodes",
        bn = x,
        check = as.integer(leaf),
        PACKAGE = "bnlearn")

}#ROOT.NODES.BACKEND

# get the parents of a node.
parents.backend = function(arcs, node, undirected = FALSE) {

  if (!undirected)
    arcs[(arcs[, "to"] == node) & which.directed(arcs), "from"]
  else
    arcs[(arcs[, "to"] == node), "from"]

}#PARENTS.BACKEND

mb.fitted = function(x, node) {

  .Call("fitted_mb",
        bn = x,
        target = node,
        PACKAGE = "bnlearn")

}#MB.FITTED

# backend of nparams, the "get the number of parameters of a
# discrete bayesian network" function. If real = TRUE this
# function returns the number of _independent_ parameters
# (one parameter of each set is set by the constraint by
# the condition \sum \pi_{i} = 1).
nparams.discrete = function(x, data, real = FALSE, debug = FALSE) {

  .Call("nparams_dnet",
        graph = x,
        data = data,
        real = real,
        debug = debug,
        PACKAGE = "bnlearn")

}#NPARAMS.DISCRETE

nparams.discrete.node = function(node, x, data, real) {

  .Call("nparams_dnode",
        graph = x,
        node = node,
        data = data,
        real = real,
        PACKAGE = "bnlearn")

}#NPARAMS.DISCRETE.NODE

nparams.gaussian = function(x, debug = FALSE) {

  .Call("nparams_gnet",
        graph = x,
        debug = debug,
        PACKAGE = "bnlearn")

}#NPARAMS.GAUSSIAN

nparams.gaussian.node = function(node, x) {

  .Call("nparams_gnode",
        graph = x,
        node = node,
        PACKAGE = "bnlearn")

}#NPARAMS.GAUSSIAN.NODE

nparams.fitted = function(x, debug = FALSE) {

  .Call("fitted_nparams",
        bn = x,
        debug = debug,
        PACKAGE = "bnlearn")

}#NPARAMS.FITTED

# return the skeleton of a graph.
dag2ug.backend = function(x, moral = FALSE, debug = FALSE) {

  nodes = names(x$nodes)

  arcs = .Call("dag2ug",
               bn = x,
               moral = moral,
               debug = debug,
               PACKAGE = "bnlearn")

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
               nodes = ordering,
               PACKAGE = "bnlearn")

  # check that the new graph is still acyclic.
  if (!is.acyclic(arcs = arcs, nodes = ordering))
    stop("this complete orientation of the graph is not acyclic.")

  return(arcs)

}#PDAG2DAG.BACKEND

# reconstruct the equivalence class of a network.
cpdag.backend = function(x, debug = FALSE) {

  nodes = names(x$nodes)

  amat = .Call("cpdag",
               arcs = x$arcs,
               nodes = nodes,
               debug = debug,
               PACKAGE = "bnlearn")

  # update the arcs of the network.
  x$arcs = amat2arcs(amat, nodes)

  # update the network structure.
  x$nodes = cache.structure(nodes, amat = amat, debug = debug)

  return(x)

}#CPDAG.BACKEND

# reconstruct the arc set of the equivalence class of a network.
cpdag.arc.backend = function(nodes, arcs, debug = FALSE) {

  amat = .Call("cpdag",
               arcs = arcs,
               nodes = nodes,
               debug = debug,
               PACKAGE = "bnlearn")

  return(amat2arcs(amat, nodes))

}#CPDAG.ARCS.BACKEND

# apply random arc operators to the graph.
perturb.backend = function(network, iter, nodes, amat, whitelist,
    blacklist, debug = FALSE) {

  # use a safety copy of the graph.
  new = network
  # remember the nodes whose score has to be recomputed.
  updates = character(0)

  # use a 'for' loop instead of a 'repeat' to avoid the threat of
  # infinite loops (due to the lack of legal operations).
  for (i in seq_len(3 * iter)) {

    to.be.added = arcs.to.be.added(amat, nodes, blacklist)
    to.be.dropped = arcs.to.be.dropped(new$arcs, whitelist)
    to.be.reversed = arcs.to.be.reversed(new$arcs, blacklist)

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
        debug = debug,
        PACKAGE = "bnlearn")

}#STRUCTURAL.HAMMING.DISTANCE

# hamming distance backend.
hamming.distance = function(learned, true, debug = FALSE) {

  .Call("shd",
        learned = skeleton(learned),
        golden = skeleton(true),
        debug = debug,
        PACKAGE = "bnlearn")

}#HAMMING.DISTANCE

# backend for extracting v-structures from a network.
vstructures = function(x, arcs, debug = FALSE) {

  .Call("vstructures",
        x = x,
        arcs = arcs,
        debug = debug,
        PACKAGE = "bnlearn")

}#VSTRUCTURES

# test equality of two graphs (nodes and arcs).
equal.backend = function(target, current) {

  .Call("all_equal",
        target = target,
        current = current,
        PACKAGE = "bnlearn")

}#EQUAL.BACKEND

# blacklists based on tiers and orderings.
tiers.backend = function(nodes, debug = FALSE) {

  .Call("tiers",
        nodes = nodes,
        debug = debug,
        PACKAGE = "bnlearn")

}#TIERS.BACKEND

cpdag.extension = function(x, debug = FALSE) {

  nodes = names(x$nodes)

  # update the arcs of the network.
  x$arcs = cpdag.arc.extension(arcs = x$arcs, nodes = nodes, debug = debug)
  # update the network structure.
  x$nodes = cache.structure(nodes, arcs = x$arcs, debug = debug)

  return(x)

}#CPDAG.EXTENSION

# backend to get a DAG out of a CPDAG (still in the same equivalence class).
cpdag.arc.extension = function(arcs, nodes, debug = FALSE) {

  .Call("pdag_extension",
        arcs = arcs,
        nodes = nodes,
        debug = debug,
        PACKAGE = "bnlearn")

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

  # construct the upper closure of the query nodes.
  upper.closure = schedule(bn, start = c(x, y, z), reverse = TRUE)
  ucgraph = subgraph.backend(bn, upper.closure)
  # get its moral graph.
  mgraph = dag2ug.backend(ucgraph, moral = TRUE)
  # remove z and incident arcs to block paths.
  upper.closure = upper.closure[!(upper.closure %in% z)]
  mgraph = subgraph.backend(mgraph, upper.closure)
  # look for a path between x and y.
  connected = has.path(x, y, nodes = upper.closure,
                amat = arcs2amat(mgraph$arcs, upper.closure))

  return(!connected)

}#UPPER.CLOSURE
