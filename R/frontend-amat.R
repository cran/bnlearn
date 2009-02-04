
# build an adjacency matrix from a graph.
amat = function(x) {

  # check x's class.
  check.bn(x)

  arcs2amat(x$arcs, names(x$nodes))

}#AMAT

# rebuild the network structure using a new adjacency matrix.
"amat<-" <- function(x, debug = FALSE, value) {

  # check x's class.
  check.bn(x)
  # a node is needed.
  if (missing(value))
    stop("no adjacency matrix specified.")
  # check the adjacency matrix.
  check.amat(amat = value, nodes = names(x$nodes))

  # update the arcs of the network.
  x$arcs = amat2arcs(value, names(x$nodes))
  # update the network structure.
  x$nodes = cache.structure(names(x$nodes), 
              amat = as.integer(value), debug = debug)

  x

}#AMAT<-

