
# original code from the deal package, released under "GPLv2 or later"
# with copyright "2002  Susanne Gammelgaard Bottcher, Claus Dethlefsen".

# write the model formula of an object of class 'bn' or 'bn.fit'.
# (ported from the deal package)
formula.backend = function(x) {

  if (is(x, "bn"))
    base = x$nodes
  else
    base = x

  paste(sapply(topological.ordering(x),
    function(node) {

      paste("[", node, ifelse(length(base[[node]]$parents) > 0, "|", ""),
        paste(base[[node]]$parents, sep="", collapse=":"), "]", sep = "")

    } ), collapse = "")

}#FORMULA.BACKEND

# create an object of class 'bn' from a model formula.
# (ported from the deal package)
model2network.backend = function(modelstring, node.order = NULL, debug = FALSE) {

  if (debug)
    cat("* processing model string:\n ", modelstring, "\n")

  # split the model string into the individual node strings.
  # the first token of the first split is always empty because the model
  # string begins with a "["; remove it.
  tokens = strsplit(strsplit(modelstring,"\\[")[[1]][-1],"\\]")

  node.set = character(0)
  arcs = vector(length(tokens), mode = "list")

  for (i in seq_along(tokens)) {

    if (debug)
      cat("  > processing node string:", tokens[[i]], "\n")

    # separate the node from its parents.
    nodes =  strsplit(tokens[[i]][1], "\\|")[[1]]
    # update the node set.
    node.set = c(node.set, nodes[1])

    if (length(nodes) > 1) {

      # if there are parents, separate them and create the corresponding arcs.
      parents = strsplit(nodes[2], ":")[[1]]

      if (debug)
        cat("    > found parents '", parents, "' for node", nodes[1], "\n")

      arcs[[i]] = cbind(from = parents, to = nodes[1])

    }#THEN
    else {

      # if there are no parents return an empty string, so that the returned
      # list can be coerced to a matrix even for an empty network.
      arcs[[i]] = cbind(from = character(0), to = character(0))

    }#THEN

  }#FOR

  # create an empty network structure.
  if (is.null(node.order)) {

    res = empty.graph.backend(node.set)

  }#THEN
  else {

    # check that the new network contains the same nodes as the node ordering.
    if (!setequal(node.set, node.order))
      stop("the nodes in the node ordering are different from those in the network.")

    res = empty.graph.backend(node.order)

  }#ELSE

  # collate and update the arcs of the network.
  res$arcs = do.call(rbind, arcs)
  # then check the arcs.
  check.arcs(res$arcs, nodes = node.set)
  # check whether the network is acyclic.
  if (!is.acyclic(nodes = node.set, arcs = res$arcs, debug = debug))
    stop("the specified network contains cycles.")
  # check whether the network is completely directed.
  if (is.pdag(arcs = res$arcs, nodes = node.set))
    stop("the graph is only partially directed.")

  # update the precomputed information on the network structure.
  if (is.null(node.order))
    res$nodes = cache.structure(sort(node.set), arcs = res$arcs, debug = debug)
  else
    res$nodes = cache.structure(node.order, arcs = res$arcs, debug = debug)

  return(res)

}#MODEL2NETWORK.BACKEND

