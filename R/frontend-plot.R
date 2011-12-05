
# generic plot of an object of class 'bn' using graphivz.
graphviz.plot = function(x, highlight = NULL, layout = "dot", shape = "circle",
    main = NULL, sub = NULL) {

  # check x's class.
  check.bn.or.fit(x)

  if (is(x, "bn")) {

    nodes = names(x$nodes)
    arcs = x$arcs

  }#THEN
  else {

    nodes = names(x)
    arcs = fit2arcs(x)

  }#ELSE

  graphviz.backend(nodes = nodes, arcs = arcs, highlight = highlight,
    layout = layout, shape = shape, main = main, sub = sub)

}#GRAPHVIZ.PLOT

# plot a graph with arcs formatted according to their own strength.
strength.plot = function(x, strength, threshold, cutpoints, highlight = NULL,
    layout = "dot", shape = "circle", main = NULL, sub = NULL, debug = FALSE) {

  # check x's class.
  check.bn(x)
  # check the strength parameter.
  check.bn.strength(strength, nodes = names(x$nodes))
  # check the strength threshold.
  threshold = check.threshold(threshold, strength)

  # check and match the strength coefficients.
  str = match.arcs.and.strengths(arcs = x$arcs, nodes = names(x$nodes),
          strengths = strength)

  # compute arc weights from the bn.strength object.
  arc.weights = strength2lwd(strength = str, threshold = threshold,
                  cutpoints = cutpoints, mode = attr(strength, "mode"),
                  arcs = x$arcs, debug = debug)

  graphviz.backend(nodes = names(x$nodes), arcs = x$arcs,
    highlight = highlight, arc.weights = arc.weights,
    layout = layout, shape = shape, main = main, sub = sub)

}#STRENGTH.PLOT

# plot method for class 'bn'.
plot.bn = function(x, ylim = c(0, 600), xlim = ylim, radius = 250, arrow = 35,
  highlight = NULL, color = "red", ...) {

  meaningless.arguments = c("xlab", "ylab", "axes")

  if (any(names(alist(...))) %in% meaningless.arguments)
    warning("arguments ", paste(meaningless.arguments, collapse = ", "),
      "will be silently ignored.")

  # check x's class.
  check.bn(x)

  if (!is.null(highlight)) {

    if (!all(highlight %in% names(x$nodes)))
      stop("invalid node label.")

    if (!(color %in% colors()) && !is.numeric(color))
      stop("invalid highlight color.")

  }#THEN
  else {

    # if there is nothing to highlight, everything is black.
    col = "black"

  }#ELSE

  # match the arguments and drop meaningless ones
  dots = eval(list(...))

  useless.ones = c('xlab', 'ylab', 'axes', 'type', 'col')
  sapply(names(dots),
    function(name) {

      if (name %in% useless.ones)
        warning("overriding the ", name, " parameter.",
          call. = FALSE)

    })

  dots[['xlab']] = dots[['ylab']] = ""
  dots[['axes']] = FALSE
  dots[['type']] = "n"
  dots[['col']] = "black"

  # create an empty plot.
  # par(oma = rep(0, 4), mar = rep(0, 4), mai = rep(0, 4),
  #   plt = c(0.06, 0.94, 0.12, 0.88))

  do.call(plot, c(list(x = 0, y = 0, xlim = xlim, ylim = ylim), dots))

  # set the stepping (in radiants) and the center of the plot.
  unit = 2 * pi / length(x$nodes)
  xc = mean(xlim)
  yc = mean(ylim)

  # compute the coordinates of the nodes
  coords = matrix(c(xc + radius * cos(1:length(x$nodes) * unit + pi/2),
                  yc + radius * sin(1:length(x$nodes) * unit + pi/2)),
                  dimnames = list(names(x$nodes), c("x" , "y")),
                  ncol = 2, byrow = FALSE)

  # draw arcs.

  if (nrow(x$arcs) > 0) apply (x$arcs, 1,
    function(a) {

      y = (coords[a[2],] - coords[a[1],]) *
            (1 - arrow / sqrt(sum((coords[a[2],] - coords[a[1],])^2))) +
            coords[a[1],]

      # if there's something to highlight, set the color according to
      # the nature of the "highlight" parameter.
      if (!is.null(highlight)) {

        if ((any(a %in% highlight) && (length(highlight) == 1)) ||
             (all(a %in% highlight) && (length(highlight) > 1)))
          col = color
        else
          col = "black"

      }#THEN

      if (is.listed(x$arcs, a[c(2, 1)]))
        length = 0
      else
        length = 0.20

      arrows(signif(coords[a[1], "x"]), signif(coords[a[1], "y"]),
             signif(y[1]), signif(y[2]), angle = 15, length = length,
             col = col)

    })

  # draw the nodes of the graph.
  points(coords, pch = 21, cex = 8, bg = "white")

  # add the names of the nodes
  for (i in 1:length(x$nodes)) {

    # if there's something to highlight, set the color according to
    # the nature of the "highlight" paramter.
    if (!is.null(highlight)) {

      if (names(x$nodes)[i] %in% highlight)
        underlined(coords[i, 1], coords[i, 2], names(x$nodes)[i], col = color)
      else
        text(coords[i, 1], coords[i, 2], names(x$nodes)[i], col = "black")

    }#THEN
    else {

      text(coords[i, 1], coords[i, 2], names(x$nodes)[i], col = "black")

    }#ELSE

  }#FOR

}#PLOT.BN


