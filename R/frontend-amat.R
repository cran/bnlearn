
# build an adjacency matrix from a graph.
amat = function(x) {

  # check x's class.
  check.bn.or.fit(x)

  if (is(x, "bn"))
    arcs2amat(x$arcs, names(x$nodes))
  else
    arcs2amat(fit2arcs(x), names(x))

}#AMAT

# rebuild the network structure using a new adjacency matrix.
"amat<-" = function(x, ignore.cycles = FALSE, debug = FALSE, value) {

  # check x's class.
  check.bn(x)
  # a node is needed.
  if (missing(value))
    stop("no adjacency matrix specified.")
  # check the adjacency matrix.
  value = check.amat(amat = value, nodes = names(x$nodes))

  # update the arcs of the network.
  x$arcs = amat2arcs(value, names(x$nodes))

  # check whether the the graph is acyclic.
  if (!ignore.cycles)
    if (!is.acyclic(nodes = names(x$nodes), arcs = x$arcs, debug = debug))
      stop("the specified network contains cycles.")

  # update the network structure.
  x$nodes = cache.structure(names(x$nodes),
              amat = as.integer(value), debug = debug)

  return(x)

}#AMAT<-

