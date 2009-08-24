
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

# get the root nodes of a network.
rootnodes.backend = function(arcs, nodes) {

  nodes[!(nodes %in% unique(arcs[, "to"]))]

}#ROOTNODES.BACKEND

# get the leaf nodes of a network.
leafnodes.backend = function(arcs, nodes) {

  nodes[!(nodes %in% unique(arcs[, "from"]))]

}#LEAFNODES.BACKEND

# get the parents of a node.
parents.backend = function(arcs, node, undirected = FALSE) {

  if (!undirected)
    arcs[(arcs[, "to"] == node) & which.directed(arcs), "from"]
  else
    arcs[(arcs[, "to"] == node), "from"]

}#PARENTS.BACKEND

# get the children of a node.
children.backend = function(arcs, node, undirected = FALSE) {

  if (!undirected)
    arcs[(arcs[, "from"] == node) & which.directed(arcs), "to"]
  else
    arcs[(arcs[, "from"] == node), "to"]

}#CHILDREN.BACKEND

# get the markov blanket of a node.
mb.backend = function(arcs, node) {

  mb = c(nbr.backend(arcs, node),
      unlist(sapply(children.backend(arcs, node),
        function(child) {

          parents.backend(arcs, node)

        }), use.names = FALSE))

  unique(mb[mb != node])

}#MB.BACKEND

# backend of nparams, the "get the number of parameters of a
# discrete bayesian network" function. If real = TRUE this
# function returns the number of _independent_ parameters
# (one parameter of each set is set by the constraint by
# the condition \sum \pi_{i} = 1).
nparams.discrete = function(x, data, real = FALSE) {

  sapply(nodes(x), nparams.discrete.node, x = x, data = data, real = real)

}#NPARAMS.DISCRETE

nparams.discrete.node = function(node, x, data, real) {

  .Call("nparams_dnode",
        graph = x,
        node = node,
        data = data,
        real = as.integer(real),
        PACKAGE = "bnlearn")

}#NPARAMS.DISCRETE.NODE

nparams.gaussian = function(x) {

  sapply(nodes(x), nparams.gaussian.node, x = x)

}#NPARAMS.GAUSSIAN

nparams.gaussian.node = function(node, x) {

  .Call("nparams_gnode",
        graph = x,
        node = node,
        PACKAGE = "bnlearn")

}#NPARAMS.GAUSSIAN.NODE

# backend for neighbourhood detection.
nbr.backend = function(arcs, node) {

  # this includes neighbours with undirected arcs.
  unique(c(arcs[arcs[, "from"] == node, "to"], arcs[arcs[, "to"] == node, "from"]))

}#NBR.BACKEND

# apply random arc operators to the graph.
perturb.backend = function(network, iter, nodes, amat, whitelist,
    blacklist, debug = FALSE) {

  # use a safety copy of the graph.
  new = network
  # remember the nodes whose score has to be recomputed.
  updates = character(0)
  # rebuild the blacklist as an adjacency matrix.
  if (!is.null(blacklist))
    blmat = arcs2amat(blacklist, nodes)
  else
    blmat = NULL

  # use a 'for' loop instead of a 'repeat' to avoid the threat of
  # infinite loops (due to the lack of legal operations).
  for (i in seq_len(3 * iter)) {

    to.be.added = arcs.to.be.added(amat, nodes, blmat)
    to.be.dropped = arcs.to.be.dropped(new$arcs, whitelist)
    to.be.reversed = arcs.to.be.reversed(new$arcs, blacklist)

    # no more operation to do.
    if (iter == 0) break;

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

