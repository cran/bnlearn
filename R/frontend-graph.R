
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
acyclic = function(x, directed, debug = FALSE) {

  # check x's class.
  check.bn.or.fit(x)
  # check debug.
  check.logical(debug)

  # fitted bayesian networks are always acylic.
  if (class(x) == "bn.fit")
    return(TRUE)

  if (missing(directed)) {

    is.acyclic(x$arcs, names(x$nodes), debug = debug)

  }#THEN
  else {

    # check directed.
    check.logical(directed)

    is.acyclic.backend(arcs = x$arcs, nodes = names(x$nodes),
      directed = directed, debug = debug)

  }#ELSE

}#ACYCLIC

# check if a graph is directed.
directed = function(x) {

  # check x's class.
  check.bn.or.fit(x)

  # fitted bayesian networks are always directed.
  if (class(x) == "bn.fit")
    return(TRUE)

  is.dag(x$arcs, names(x$nodes))

}#DIRECTED

# check if there's a path between two specific nodes.
path = function(x, from, to, direct = TRUE, underlying.graph = FALSE,
    debug = FALSE) {

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

  if (class(x) == "bn") {

    nodes = names(x$nodes)
    amat = arcs2amat(x$arcs, nodes)

  }#THEN
  else {

    nodes = names(x)
    amat = arcs2amat(fit2arcs(x), nodes)

  }#ELSE

  has.path(from, to, nodes = nodes, amat = amat, exclude.direct = !direct,
    underlying.graph = underlying.graph, debug = debug)

}#PATH

# return the partial node ordering implied by the graph structure.
node.ordering = function(x, debug = FALSE) {

  # check x's class.
  check.bn.or.fit(x)
  # check debug.
  check.logical(debug)
  # no model string if the graph is partially directed.
  if (class(x) == "bn")
    if (is.pdag(x$arcs, names(x$nodes)))
      stop("the graph is only partially directed.")

  schedule(x, debug = debug)

}#NODE.ORDERING

# generate a valid blacklist from a partial node ordering.
ordering2blacklist = function(nodes) {

  if (class(nodes) %in% c("bn", "bn.fit")) {

    nodes = schedule(nodes)

  }#THEN

  # check the node labels.
  check.nodes(nodes, min.nodes = 3)

  o2b.backend(nodes)

}#ORDERING2BLACKLIST

# return the skeleton of a (partially) directed graph
skeleton = function(x) {

  # check x's class.
  check.bn(x)

  dag2ug.backend(x, moral = FALSE)

}#SKELETON

# return a complete orientation of a graph.
pdag2dag = function(x, ordering) {

  # check x's class.
  check.bn(x)
  # check the node ordering.
  check.nodes(ordering, graph = x, min.nodes = length(x$nodes),
    max.nodes = length(x$nodes))

  pdag2dag.backend(x, ordering)

}#PDAG2DAG

# test the equality of two network structures.
all.equal.bn = function(target, current, ...) {

  # check the class of target and current.
  check.bn(target)
  check.bn(current)
  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  equal.backend(target, current)

}#ALL.EQUAL.BN

# compare two bayesian network structures.
compare = function(target, current, debug = FALSE) {

  result = TRUE

  # check both objects' class.
  check.bn(target)
  check.bn(current)
  # check debug.
  check.logical(debug)

  compare.backend(target = target, current = current, debug = debug)

}#COMPARE

