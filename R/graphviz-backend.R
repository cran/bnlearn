
# unified backend for the graphviz calls.
graphviz.backend = function(nodes, arcs, highlight = NULL, groups,
    arc.weights = NULL, layout = "dot", shape = "circle", main = NULL,
    sub = NULL, render = TRUE) {

  node.shapes = c("ellipse", "circle", "rectangle")
  highlight.params = c("nodes", "arcs", "col", "fill", "lwd", "lty", "textCol")
  highlighting = FALSE

  # sanitize the layout (to be passed to layoutGraph).
  check.label(layout, choices = graphviz.layouts, argname = "graph layout")
  # sanitize nodes' shape.
  if (!is.string(shape))
    stop("node shape must be a character string.")
  if (shape %!in% node.shapes)
    stop("valid node schapes are:", paste0(" '", node.shapes, "'"), ".")
  # sanitize arc weights.
  if (!is.null(arc.weights)) {

    if (!is.numeric(arc.weights))
      stop("arc weights must be numeric values.")
    if (length(arc.weights) != nrow(arcs))
      stop("mismatch between the number of weights and the number of arcs.")

  }#THEN
  # sanitize the highlighting parameters (to be saved in *RenderInfo()).
  if (!is.null(highlight) || length(highlight) > 0) {

    highlighting = TRUE

    if (!is.list(highlight) || any(names(highlight) %!in% highlight.params))
      stop("highlight must be a list with a subset of the following",
           " elements:", paste0(" '", highlight.params, "'"), ".")

    if ("nodes" %in% names(highlight))
      check.nodes(highlight$nodes, graph = nodes)

    if ("arcs" %in% names(highlight)) {

      highlight$arcs = check.arcs(highlight$arcs, nodes = nodes)

      # disregard empty arc sets.
      if (nrow(highlight$arcs) == 0)
        highlight$arcs = NULL

    }#THEN

    if ("col" %in% names(highlight))
      check.colour(highlight$col)
    else
      highlight$col = "red"

    if ("fill" %in% names(highlight)) {

      if ("nodes" %!in% names(highlight))
        warning("no node to apply the 'fill' color to, ignoring.")

      check.colour(highlight$fill)

    }#THEN
    else
      highlight$fill = "transparent"

    if ("lwd" %in% names(highlight)) {

      if ("arcs" %!in% names(highlight))
        warning("no arc to apply the 'lwd' setting to, ignoring.")

      if (!is.positive(highlight$lwd))
        stop("the line width must be a positive number.")

    }#THEN

    if ("lty" %in% names(highlight)) {

      if ("arcs" %!in% names(highlight))
        warning("no arc to apply the 'lty' setting to, ignoring.")

      check.lty(highlight$lty)

    }#THEN

    if ("textCol" %in% names(highlight)) {

      if ("nodes" %!in% names(highlight))
        warning("no node to apply the 'textColor' color to, ignoring.")

      check.colour(highlight$textCol)

    }#THEN
    else
      highlight$textCol = "black"

  }#THEN

  # create the graphNEL object from the bn object.
  graph.obj = new("graphNEL", nodes = nodes, edgeL = arcs2elist(arcs, nodes),
                edgemode = 'directed')

  # set the title and the subtitle.
  graph::graphRenderInfo(graph.obj)[["main"]] = main
  graph::graphRenderInfo(graph.obj)[["sub"]] = sub

  # set graph layout and global parameters.
  if (shape %in% c("ellipse", "rectangle")) {

    # nodes have elliptic/rectangular shapes whose width is determined by the
    # length of the respective labels.
    attrs = list(node = list(fixedsize = FALSE))
    node.attrs = list(shape = rep(shape, length(nodes)))
    names(node.attrs$shape) = nodes

  }#THEN
  else if (shape == "circle") {

    # nodes have circular shapes.
    attrs = node.attrs = list()

  }#THEN

  # create subgraphs from the node groups.
  if (!missing(groups) && !is.null(groups)) {

    # check the node groups.
    check.node.groups(groups, graph = nodes)

    subGList = lapply(groups, function(g) {
      list(graph = graph::subGraph(g, graph.obj), cluster = TRUE)
    })

  }#THEN
  else {

    subGList = NULL

  }#ELSE

  graph.plot = Rgraphviz::layoutGraph(graph.obj, subGList = subGList,
                 attrs = attrs, nodeAttrs = node.attrs, layoutType = layout)

  # default arcs to solid lines of width 1; otherwise, when we set it for some
  # arcs via weights or highlight the remaining arcs have it set to NA.
  graph::edgeRenderInfo(graph.plot)[["lty"]] = "solid"
  graph::edgeRenderInfo(graph.plot)[["lwd"]] = 1

  # kill the arrowheads of undirected arcs.
  u = names(which(graph::edgeRenderInfo(graph.plot)[['direction']] == "both"))

  graph::edgeRenderInfo(graph.plot)[["arrowhead"]][u] = "none"
  graph::edgeRenderInfo(graph.plot)[["arrowtail"]][u] = "none"

  # change arc line width according to arc weights.
  if (!is.null(arc.weights) && (length(arc.weights) > 0)) {

    to.weight = apply(arcs, 1, paste, collapse = "~")

    for (i in 1:length(to.weight)) {

      # plot an arc as a dotted line if it has a negative weight
      # (i.e. it's removal would improve the goodness of fit).
      if (arc.weights[i] > 0)
        graph::edgeRenderInfo(graph.plot)[["lwd"]][to.weight[i]] = arc.weights[i]
      else
        graph::edgeRenderInfo(graph.plot)[["lty"]][to.weight[i]] = "dashed"

    }#THEN

  }#THEN

  # apply the requested highlighting.
  if (highlighting) {

    if ("nodes" %in% names(highlight)) {

      graph::nodeRenderInfo(graph.plot)[["col"]][highlight$nodes] = highlight$col
      graph::nodeRenderInfo(graph.plot)[["fill"]][highlight$nodes] = highlight$fill
      graph::nodeRenderInfo(graph.plot)[["textCol"]][highlight$nodes] = highlight$textCol

    }#THEN

    if ("arcs" %in% names(highlight)) {

      to.highlight = apply(highlight$arcs, 1, paste, collapse = "~")

      graph::edgeRenderInfo(graph.plot)[["col"]][to.highlight] = highlight$col

      if ("lwd" %in% names(highlight))
        graph::edgeRenderInfo(graph.plot)[["lwd"]][to.highlight] = highlight$lwd

      # note that this overrides the changes made according to arc weights.
      if ("lty" %in% names(highlight))
        graph::edgeRenderInfo(graph.plot)[["lty"]][to.highlight] = highlight$lty

    }#THEN

  }#THEN

  if (render) {

    # do the actual plotting.
    if (nrow(arcs) > 0) {

      Rgraphviz::renderGraph(graph.plot)

    }#THEN
    else {

      # use a NOP function to render arcs, the default renderEdges() raises an
      # error if the graph is empty.
      Rgraphviz::renderGraph(graph.plot, drawEdges = function(x) {})

    }#ELSE

  }#THEN

  # return (invisibly) the graph object, to allow further customizations.
  invisible(graph.plot)

}#GRAPHVIZ.BACKEND

