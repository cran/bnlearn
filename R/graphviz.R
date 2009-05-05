
# unified backend for the graphviz calls.
graphviz.backend = function(nodes, arcs, highlight = NULL, arc.weights = NULL,
    layout = "dot", main = NULL, sub = NULL) {

  graphviz.layouts = c("dot", "neato", "twopi", "circo", "fdp")
  highlight.params = c("nodes", "arcs", "col", "fill")
  highlighting = FALSE

  # check whether graphviz is loaded.
  if (!("Rgraphviz" %in% loadedNamespaces()))
    stop("this function requires Rgraphviz.")

  # sanitize the layout (to be passed to layoutGraph).
  if (!(layout %in% graphviz.layouts))
    stop(paste(c("valid layout schemes are:", graphviz.layouts), collapse = " "))
  # sanitize the highlighting parameters (to be saved in *RenderInfo()).
  if (!is.null(highlight) || length(highlight) > 0) {

    highlighting = TRUE

    if (!is.list(highlight) || !all(names(highlight) %in% highlight.params))
      stop(paste(c("highlight must be a list with at least one of the",
             "following elements:", highlight.params), collapse = " "))

    fake.graph = list(nodes = structure(nodes, names = nodes))

    if ("nodes" %in% names(highlight))
      check.nodes(highlight$nodes, graph = fake.graph)

    if ("arcs" %in% names(highlight))
      highlight$arcs = check.arcs(highlight$arcs, graph = fake.graph)

    if ("col" %in% names(highlight))
      check.colour(highlight$col)
    else
      highlight$col = "red"

    if ("fill" %in% names(highlight)) {

      if (!("nodes" %in% names(highlight)))
        warning("no node to apply the 'fill' color to, ignoring.")

      check.colour(highlight$fill)

    }#THEN
    else
      highlight$fill = "transparent"

  }#THEN

  # create the graphAM object from the bn object.
  graph.obj = new("graphAM", adjMat = arcs2amat(arcs, nodes),
    edgemode = 'directed')

  # dump the global graphical settings.
  graph.par.dump = graph.par()

  # set the title and the subtitle.
  graph.par(list(graph = list(main = main, sub = sub)))

  # set graph layout and global parameters.
  graph.plot = layoutGraph(graph.obj, layoutType = layout)

  # kill the arroheads of undirected arcs.
  u = names(which(edgeRenderInfo(graph.plot)[['direction']] == "both"))

  edgeRenderInfo(graph.plot)[["arrowhead"]][u] = "none"
  edgeRenderInfo(graph.plot)[["arrowtail"]][u] = "none"

  # apply the requested highlighting.
  if (highlighting) {

    if ("nodes" %in% names(highlight)) {

      nodeRenderInfo(graph.plot)[["col"]][highlight$nodes] = highlight$col
      nodeRenderInfo(graph.plot)[["fill"]][highlight$nodes] = highlight$fill

    }#THEN

    if ("arcs" %in% names(highlight)) {

      to.highlight = apply(highlight$arcs, 1, paste, collapse = "~")

      edgeRenderInfo(graph.plot)[["col"]][to.highlight] = highlight$col

    }#THEN

  }#THEN

  # change arc line width according to arc weights.
  if (!is.null(arc.weights)) {

    to.weight = apply(arcs, 1, paste, collapse = "~")

    for (i in 1:length(to.weight)) {

      # plot an arc as a dotted line if it has a negative weight
      # (i.e. it's removal would improve the goodness of fit).
      if (arc.weights[i] > 0)
        edgeRenderInfo(graph.plot)[["lwd"]][to.weight[i]] = arc.weights[i]
      else
        edgeRenderInfo(graph.plot)[["lty"]][to.weight[i]] = "dashed"

    }#THEN

  }#THEN

  # do the actual plotting
  renderGraph(graph.plot)

  # restore the original global graphical settings.
  graph.par(graph.par.dump)

  # return the graph object, to allow further customizations.
  return(graph.plot)

}#GRAPHVIZ.BACKEND

