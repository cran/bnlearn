
graphviz.plot = function(x, highlight = NULL, layout = "dot") {

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

    if ("nodes" %in% names(highlight))
      check.nodes(highlight$nodes, graph = x)

    if ("arcs" %in% names(highlight))
      highlight$arcs = check.arcs(highlight$arcs, graph = x)

    if ("col" %in% names(highlight))
      col2rgb(highlight$col)
    else
      highlight$col = "red"

    if ("fill" %in% names(highlight))
      col2rgb(highlight$col)
    else
      highlight$fill = "transparent"

  }#THEN

  # create the graphAM object from the bn object.
  graph.obj = new("graphAM", adjMat = amat(x), edgemode = 'directed')

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

      arcs = apply(highlight$arcs, 1, paste, collapse = "~")

      edgeRenderInfo(graph.plot)[["col"]][arcs] = highlight$col

    }#THEN

  }#THEN

  # do the actual plotting
  renderGraph(graph.plot)

  # return the graph object, to allow further customizations.
  return(graph.plot)

}#GRAPHVIZ.PLOT

