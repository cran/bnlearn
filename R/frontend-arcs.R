
# return the arcs in the graph.
arcs = function(x) {

  # check x's class.
  check.bn(x)

  x$arcs

}#ARCS

# rebuild the network structure using a new set fo arcs.
"arcs<-" <- function(x, ignore.cycles = FALSE, debug = FALSE, value) {

  # check x's class.
  check.bn(x)
  # a set of arcs is needed.
  if (missing(value))
    stop("no arc specified.")
  # sanitize the set of arcs.
  value = check.arcs(value, graph = x)
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
  check.bn(x)

  x$arcs[!which.undirected(x$arcs), , drop = FALSE]

}#DIRECTED.ARCS

# return the undirected arcs in the graph.
undirected.arcs = function(x) {

  # check x's class.
  check.bn(x)

  x$arcs[which.undirected(x$arcs), , drop = FALSE]

}#UNDIRECTED.ARCS

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

