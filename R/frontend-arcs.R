
# return the arcs in the graph.
arcs = function(x) {

  # check x's class.
  check.bn.or.fit(x)

  if (is(x, "bn"))
    x$arcs
  else
    fit2arcs(x)

}#ARCS

# rebuild the network structure using a new set fo arcs.
"arcs<-" = function(x, ignore.cycles = FALSE, debug = FALSE, value) {

  # check x's class.
  check.bn(x)
  # a set of arcs is needed.
  if (missing(value))
    stop("no arc specified.")
  # sanitize the set of arcs.
  value = check.arcs(value, nodes = names(x$nodes))
  # check whether the the graph is acyclic.
  if (!ignore.cycles)
    if (!is.acyclic(nodes = names(x$nodes), arcs = value, debug = debug))
      stop("the specified network contains cycles.")

  # update the arcs of the network.
  x$arcs = value
  # update the network structure.
  x$nodes = cache.structure(names(x$nodes), arcs = x$arcs, debug = debug)

  return(x)

}#ARCS<-

# return the directed arcs in the graph.
directed.arcs = function(x) {

  # check x's class.
  check.bn.or.fit(x)

  if (is(x, "bn"))
    x$arcs[which.directed(x$arcs, names(x$nodes)), , drop = FALSE]
  else
    fit2arcs(x)

}#DIRECTED.ARCS

# return the undirected arcs in the graph.
undirected.arcs = function(x) {

  # check x's class.
  check.bn.or.fit(x)

  if (is(x, "bn"))
    x$arcs[which.undirected(x$arcs, names(x$nodes)), , drop = FALSE]
  else
    matrix(character(0), nrow = 0, ncol = 2,
      dimnames = list(NULL, c("from", "to")))

}#UNDIRECTED.ARCS

# return the arcs pointing to a particular node.
incoming.arcs = function(x, node) {

  # check x's class.
  check.bn.or.fit(x)
  # a valid node is needed.
  check.nodes(nodes = node, graph = x, max.nodes = 1)

  arcs = directed.arcs(x)

  arcs[arcs[, "to"] == node, , drop = FALSE]

}#INCOMING.ARCS

# return the arcs originating from a particular node.
outgoing.arcs = function(x, node) {

  # check x's class.
  check.bn.or.fit(x)
  # a valid node is needed.
  check.nodes(nodes = node, graph = x, max.nodes = 1)

  arcs = directed.arcs(x)

  arcs[arcs[, "from"] == node, , drop = FALSE]

}#OUTGOING.ARCS

# return the arcs incident on a particular node.
incident.arcs = function(x, node) {

  # check x's class.
  check.bn.or.fit(x)
  # a valid node is needed.
  check.nodes(nodes = node, graph = x, max.nodes = 1)

  arcs = arcs(x)

  arcs[(arcs[, "from"] == node) | (arcs[, "to"] == node), , drop = FALSE]

}#INCIDENT.ARCS

# return the number of arcs in the graph.
narcs = function(x) {

  # check x's class.
  check.bn.or.fit(x)

  narcs.backend(x)

}#NARCS

# set an arc direction manually.
set.arc = function(x, from, to, check.cycles = TRUE, debug = FALSE) {

  arc.operations(x = x, from = from, to = to, op = "set",
    check.cycles = check.cycles, debug = debug)

}#SET.ARC

# remove an arc from the graph.
drop.arc = function(x, from, to, debug = FALSE) {

  arc.operations(x = x, from = from, to = to, op = "drop",
    check.cycles = FALSE, debug = debug)

}#DROP.ARC

# reverse an arc in the graph.
reverse.arc = function(x, from, to, check.cycles = TRUE, debug = FALSE) {

  arc.operations(x = x, from = from, to = to, op = "reverse",
    check.cycles = check.cycles, debug = debug)

}#REVERSE.ARC

