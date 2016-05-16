
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
"amat<-" = function(x, check.cycles = TRUE, check.illegal = TRUE, debug = FALSE,
    value) {

  # check x's class.
  check.bn(x)
  # a node is needed.
  if (missing(value))
    stop("no adjacency matrix specified.")
  # check logical arguments.
  check.logical(check.cycles)
  check.logical(check.illegal)
  check.logical(debug)
  # check the adjacency matrix.
  value = check.amat(amat = value, nodes = names(x$nodes))

  # update the arcs of the network.
  x$arcs = amat2arcs(value, names(x$nodes))

  # check whether the the graph contains directed cycles.
  if (check.cycles)
    if (!is.acyclic(nodes = names(x$nodes), arcs = x$arcs, debug = debug,
           directed = TRUE))
      stop("the specified network contains cycles.")
  # check whether any arc is illegal.
  if (check.illegal) {

    illegal = which.listed(x$arcs, x$learning$illegal)

    if (any(illegal)) {

      illegal = apply(x$arcs[illegal, , drop = FALSE], 1,
                  function(x) { paste(" (", x[1], ", ", x[2], ")", sep = "")  })

      stop("the following arcs are not valid due to the parametric assumptions of the network:",
        illegal, ".")

    }#THEN

  }#THEN

  # update the network structure.
  x$nodes = cache.structure(names(x$nodes),
              amat = as.integer(value), debug = debug)

  return(x)

}#AMAT<-

