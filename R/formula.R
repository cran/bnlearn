
# original code from the deal package, released under "GPLv2 or later"
# with copyright "2002  Susanne Gammelgaard Bottcher, Claus Dethlefsen".

# write the model formula of an object of class 'bn' or 'bn.fit'.
# (ported from the deal package)
formula.backend = function(x) {

  if (is(x, "bn"))
    base = x$nodes
  else
    base = x

  paste(sapply(schedule(x),
    function(node) {

      paste("[", node, ifelse(length(base[[node]]$parents) > 0, "|", ""),
        paste(base[[node]]$parents, sep="", collapse=":"), "]", sep = "")

    } ), collapse = "")

}#FORMULA.BACKEND

# create an object of class 'bn' from a model formula.
# (ported from the deal package)
model2network.backend = function(modelstring, node.order = NULL, debug = FALSE) {

    nodes = c()

    if (debug)
      cat("* processing model string:\n ", modelstring, "\n")

    # split the model strings into the single nodes' strings.
    # the first entry of the first split in always empty because the model
    # string begins with a "["; remove it.
    st = strsplit(strsplit(modelstring,"\\[")[[1]][-1],"\\]")

    arcs = lapply(st, function(xp) {

      if (debug)
       cat("  > processing node string:", xp, "\n")

      # separate the node from its parents.
      pa =  strsplit(xp[1], "\\|")[[1]]

      assign('nodes', c(nodes, pa[1]), envir = sys.frame(-2))

      # if there are parents ...
      if (length(pa) > 1) {

        prnts = strsplit(pa[2], ":")[[1]]

        if (debug)
          cat("    > detected parents '", prnts, "' for node", pa[1], "\n")

        return(cbind(from = prnts, to = pa[1]))

      }#THEN
      else {

        # if there are no parents return an empty string, so that the returned
        # list can be coerced to a matrix even for an empty graph.
        return(cbind(from = character(0), to = character(0)))

      }#THEN

    })

    # create an empty network structure.
    if (is.null(node.order)) {

      res = empty.graph.backend(nodes)

    }#THEN
    else {

      # check that the new network contains the same nodes as the old one.
      if (!setequal(nodes, node.order))
        stop("the model string describes a different network (the nodes are different).")

      res = empty.graph.backend(node.order)

    }#ELSE

    # update the arcs of the network.
    res$arcs = do.call(rbind, arcs)
    # then check the arcs.
    check.arcs(res$arcs, nodes = nodes)
    # check whether the the graph is acyclic.
    if (!is.acyclic(nodes = nodes, arcs = res$arcs, debug = debug))
      stop("the specified network contains cycles.")

    # update the network structure.
    if (is.null(node.order)) {

      res$nodes = cache.structure(sort(nodes), arcs = res$arcs, debug = debug)

    }#THEN
    else {

      res$nodes = cache.structure(node.order, arcs = res$arcs, debug = debug)

    }#ELSE

    return(res)

}#MODEL2NETWORK

