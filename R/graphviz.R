
# unified backend for the graphviz calls.
graphviz.backend = function(nodes, arcs, highlight = NULL, arc.weights = NULL,
    layout = "dot", shape = "circle", main = NULL, sub = NULL) {

  graphviz.layouts = c("dot", "neato", "twopi", "circo", "fdp")
  node.shapes = c("ellipse", "circle")
  highlight.params = c("nodes", "arcs", "col", "fill", "lwd", "lty")
  highlighting = FALSE

  # check whether graphviz is loaded.
  if (!("Rgraphviz" %in% loadedNamespaces()))
    stop("this function requires Rgraphviz.")

  # sanitize the layout (to be passed to layoutGraph).
  if (!is.string(layout))
    stop("graph layout must be a character string.")
  if (!(layout %in% graphviz.layouts))
    stop(paste(c("valid layout schemes are:", graphviz.layouts), collapse = " "))
  # sanitize nodes' shape.
  if (!is.string(shape))
    stop("node shape must be a character string.")
  if (!(shape %in% node.shapes))
    stop(paste(c("valid node shapes are:", node.shapes), collapse = " "))
  # sanitize the highlighting parameters (to be saved in *RenderInfo()).
  if (!is.null(highlight) || length(highlight) > 0) {

    highlighting = TRUE

    if (!is.list(highlight) || !all(names(highlight) %in% highlight.params))
      stop(paste(c("highlight must be a list with at least one of the",
             "following elements:", highlight.params), collapse = " "))

    if ("nodes" %in% names(highlight))
      check.nodes(highlight$nodes, graph = nodes)

    if ("arcs" %in% names(highlight))
      highlight$arcs = check.arcs(highlight$arcs, nodes = nodes)

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

    if ("lwd" %in% names(highlight)) {

      if (!("arcs" %in% names(highlight)))
        warning("no arc to apply the 'lwd' setting to, ignoring.")

      if (!is.positive(highlight$lwd))
        stop("the line width must be a positive number.")

    }#THEN

    if ("lty" %in% names(highlight)) {

      if (!("arcs" %in% names(highlight)))
        warning("no arc to apply the 'lty' setting to, ignoring.")

      check.lty(highlight$lty)

    }#THEN

  }#THEN

  # create the graphAM object from the bn object.
  graph.obj = new("graphAM", adjMat = arcs2amat(arcs, nodes),
    edgemode = 'directed')

  # dump the global graphical settings.
  graph.par.dump = graph.par()

  # set the title and the subtitle.
  graph.par(list(graph = list(main = main, sub = sub)))

  # set graph layout and global parameters.
  if (shape == "ellipse") {

    # nodes have elliptic shapes whose width is determined by the
    # length of the respective labels.
    attrs = list(node = list(fixedsize = FALSE))
    node.attrs = list(shape = rep("ellipse", length(nodes)))
    names(node.attrs$shape) = nodes

  }#THEN
  else if (shape == "circle") {

    # nodes have circular shapes.
    attrs = node.attrs = list()

  }#THEN

  graph.plot = layoutGraph(graph.obj, attrs = attrs, nodeAttrs = node.attrs,
                 layoutType = layout)

  # kill the arroheads of undirected arcs.
  u = names(which(edgeRenderInfo(graph.plot)[['direction']] == "both"))

  edgeRenderInfo(graph.plot)[["arrowhead"]][u] = "none"
  edgeRenderInfo(graph.plot)[["arrowtail"]][u] = "none"

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

  # apply the requested highlighting.
  if (highlighting) {

    if ("nodes" %in% names(highlight)) {

      nodeRenderInfo(graph.plot)[["col"]][highlight$nodes] = highlight$col
      nodeRenderInfo(graph.plot)[["fill"]][highlight$nodes] = highlight$fill

    }#THEN

    if ("arcs" %in% names(highlight)) {

      to.highlight = apply(highlight$arcs, 1, paste, collapse = "~")

      edgeRenderInfo(graph.plot)[["col"]][to.highlight] = highlight$col

      if ("lwd" %in% names(highlight))
        edgeRenderInfo(graph.plot)[["lwd"]][to.highlight] = highlight$lwd

      # note that this overrides the changes made according to arc weights.
      if ("lty" %in% names(highlight))
        edgeRenderInfo(graph.plot)[["lty"]][to.highlight] = highlight$lty

    }#THEN

  }#THEN

  # do the actual plotting.
  if (nrow(arcs) > 0) {

    renderGraph(graph.plot)

  }#THEN
  else {

    # use a NOP function to render arcs, the default renderEdges() raises an
    # error if the graph is empty.
    renderGraph(graph.plot, drawEdges = function(x) {})

  }#ELSE

  # restore the original global graphical settings.
  graph.par(graph.par.dump)

  # return (invisibly) the graph object, to allow further customizations.
  invisible(graph.plot)

}#GRAPHVIZ.BACKEND

