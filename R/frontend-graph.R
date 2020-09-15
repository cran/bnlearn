# check if there's a path between two specific nodes.
path.exists = function(x, from, to, direct = TRUE,
    underlying.graph = FALSE, debug = FALSE) {

  # check x's class.
  check.bn.or.fit(x)
  # a valid node is needed.
  check.nodes(nodes = from, graph = x, max.nodes = 1)
  # another valid node is needed.
  check.nodes(nodes = to, graph = x, max.nodes = 1)
  # 'from' must be different from 'to'.
  if (identical(from, to))
    stop("'from' and 'to' must be different from each other.")
  # check underlying.path.
  check.logical(underlying.graph)
  # check debug.
  check.logical(debug)

  if (is(x, "bn")) {

    nodes = names(x$nodes)
    amat = arcs2amat(x$arcs, nodes)

  }#THEN
  else {

    nodes = names(x)
    amat = arcs2amat(fit2arcs(x), nodes)

  }#ELSE

  has.path(from, to, nodes = nodes, amat = amat, exclude.direct = !direct,
    underlying.graph = underlying.graph, debug = debug)

}#PATH.EXISTS

# get the number of nodes of a graph.
nnodes = function(x) {

  # check x's class.
  check.bn.or.fit(x)

  if (is(x, "bn"))
    return(length(x$nodes))
  else
    return(length(x))

}#NNODES

# get the root nodes of a graph.
root.nodes = function(x) {

  # check x's class.
  check.bn.or.fit(x)

  root.leaf.nodes(x, leaf = FALSE)

}#ROOT.NODES

# get the leaf nodes of a graph.
leaf.nodes = function(x) {

  # check x's class.
  check.bn.or.fit(x)

  root.leaf.nodes(x, leaf = TRUE)

}#LEAF.NODES

# check if a graph is acyclic.
acyclic = function(x, directed = FALSE, debug = FALSE) {

  # check x's class.
  check.bn.or.fit(x)
  # check debug and directed.
  check.logical(debug)
  check.logical(directed)

  # fitted bayesian networks are always acylic.
  if (is(x, "bn.fit"))
    return(TRUE)

  is.acyclic(arcs = x$arcs, nodes = names(x$nodes), debug = debug,
    directed = directed)

}#ACYCLIC

# check if a graph is directed.
directed = function(x) {

  # check x's class.
  check.bn.or.fit(x)

  # fitted bayesian networks are always directed.
  if (is(x, "bn.fit"))
    return(TRUE)

  is.dag(x$arcs, names(x$nodes))

}#DIRECTED

# return the partial node ordering implied by the graph structure.
node.ordering = function(x, debug = FALSE) {

  # check x's class.
  check.bn.or.fit(x)
  # check debug.
  check.logical(debug)
  # no model string if the graph is partially directed.
  if (is(x, "bn"))
    if (is.pdag(x$arcs, names(x$nodes)))
      stop("the graph is only partially directed.")

  topological.ordering(x, debug = debug)

}#NODE.ORDERING

# generate a valid blacklist from a partial node ordering.
ordering2blacklist = function(nodes) {

  if (is(nodes, "bn") || is(nodes, "bn.fit"))
    nodes = topological.ordering(nodes)

  # check the node labels.
  check.nodes(nodes)

  tiers.backend(nodes)

}#ORDERING2BLACKLIST

# generate a valid blacklist from a partial node ordering.
tiers2blacklist = function(tiers) {

  # check the node labels.
  check.node.groups(tiers)

  tiers.backend(tiers)

}#TIERS2BLACKLIST

# generate a blacklist containing all arcs between nodes in a set.
set2blacklist = function(set) {

  if (!is.string.vector(set))
    stop("'set' should be a character vector, the labels of the nodes.")

  # create the blacklist, excluding loops.
  bl = expand.grid(from = set, to = set, stringsAsFactors = FALSE)
  bl = bl[bl$from != bl$to, ]
  # reset row numbers.
  rownames(bl) = NULL

  return(as.matrix(bl))

}#SET2BLACKLIST

# return the skeleton of a (partially) directed graph.
skeleton = function(x) {

  # check x's class.
  check.bn.or.fit(x)

  # go back to the network structure if needed.
  if (is(x, "bn.fit"))
    x = bn.net(x)

  dag2ug.backend(x, moral = FALSE)

}#SKELETON

# return a complete orientation of a graph given a topological ordering.
pdag2dag = function(x, ordering) {

  # check x's class.
  check.bn(x)
  # check the node ordering.
  check.nodes(ordering, graph = x, min.nodes = length(x$nodes),
    max.nodes = length(x$nodes))

  arcs = pdag2dag.backend(x$arcs, ordering)

  # update the arcs of the network.
  x$arcs = arcs
  # update the network structure.
  x$nodes = cache.structure(nodes = names(x$nodes), arcs = arcs)

  return(x)

}#PDAG2DAG

# test the equality of two network structures.
all.equal.bn = function(target, current, ...) {

  # check the class of target and current.
  check.bn(target)
  check.bn(current)
  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  equal.backend.bn(target, current)

}#ALL.EQUAL.BN

# compare two bayesian network structures.
compare = function(target, current, arcs = FALSE) {

  # check both objects' class.
  check.bn(target)
  check.bn(current)
  # the two networks must have the same node set.
  match.bn(target, current)
  # check debug.
  check.logical(arcs)

  compare.backend(target$arcs, current$arcs,
    nodes = names(target$nodes), arcs = arcs)

}#COMPARE

# create a subgraph spanning a subset of nodes.
subgraph = function(x, nodes) {

  # check x's class.
  check.bn.or.fit(x)
  # get the network structure out of a bn.fit object.
  if (is(x, "bn.fit"))
    x = bn.net(x)
  # check the nodes of the subgraph.
  check.nodes(nodes, graph = x, max.nodes = length(x$nodes))

  # creating a new graph, so do not preserve learning information.
  subgraph.backend(x = x, nodes = nodes, preserve.learning = FALSE)

}#SUBGRAPH

