
# generic plot of an object of class 'bn' or 'bn.fit' using graphviz.
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

  # return the graph object for further customization, not NULL.
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
                  cutpoints = cutpoints, method = attr(strength, "method"),
                  arcs = x$arcs, debug = debug)

  # return the graph object for further customization, not NULL.
  graphviz.backend(nodes = names(x$nodes), arcs = x$arcs,
    highlight = highlight, arc.weights = arc.weights,
    layout = layout, shape = shape, main = main, sub = sub)

}#STRENGTH.PLOT

# plot method for class 'bn'.
plot.bn = function(x, ylim = c(0, 600), xlim = ylim, radius = 250, arrow = 35,
  highlight = NULL, color = "red", ...) {

  # check x's class.
  check.bn(x)

  if (!is.null(highlight)) {

    if (any(highlight %!in% names(x$nodes)))
      stop("invalid node label.")

    if ((color %!in% colors()) && !is.numeric(color))
      stop("invalid highlight color.")

  }#THEN
  else {

    # if there is nothing to highlight, everything is black.
    col = "black"

  }#ELSE

  # match the arguments and drop meaningless ones.
  dots = sanitize.plot.dots(eval(list(...)),
           meaningless = c("xlab", "ylab", "axes", "type", "col"))
  # set a few defaults.
  dots[['xlab']] = dots[['ylab']] = ""
  dots[['axes']] = FALSE
  dots[['type']] = "n"
  dots[['col']] = "black"

  # create an empty plot.
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

  # draw the arcs.
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

  invisible(NULL)

}#PLOT.BN

# ECDF plot for 'bn.strength' objects.
plot.bn.strength = function(x, draw.threshold = TRUE, main = NULL,
    xlab = "arc strengths", ylab = "CDF(arc strengths)", ...) {

  # check the draw.threshold flag.
  check.logical(draw.threshold)
  # check the arcs strengths.
  if (missing(x))
    stop("an object of class 'bn.strength' is required.")
  if (!is(x, "bn.strength"))
    stop("x must be a 'bn.strength' object.")
  # only arc strengths computed as confidence through bootstrap are supported.
  if (attr(x, "method") != "bootstrap")
    stop("only arc strengths obtained through bootstrap are supported.")

  # match the dots arguments.
  dots = sanitize.plot.dots(eval(list(...)), meaningless = c("xlim", "ylim"))

  # handle the title.
  if (is.null(main))
    main = paste("threshold = ", attr(x, "threshold"))

  # draw the plot.
  do.call(plot, c(list(x = ecdf(x$strength), xlim = c(-0.1, 1.1), ylim = c(0, 1),
    xlab = xlab, ylab = ylab, main = main), dots))

  if (draw.threshold)
    abline(v = attr(x, "threshold"), lty = 2)

  invisible(NULL)

}#PLOT.BN.STRENGTH

# boxplot for cross-validation losses (single run).
plot.bn.kcv = function(x, ..., main, xlab, ylab, connect = FALSE) {

  # check x's class.
  if (!is(x, "bn.kcv"))
    stop("x must be an object of class 'bn.kcv'.")

  plot.bn.kcv.list(x = structure(list(x), class = "bn.kcv.list"), ...,
    main = main, xlab = xlab, ylab = ylab, connect = connect)

  invisible(NULL)

}#PLOT.BN.KCV

# boxplot for cross-validation losses (multiple runs).
plot.bn.kcv.list = function(x, ..., main, xlab, ylab, connect = FALSE) {

  # check x's class.
  if (!is(x, "bn.kcv.list"))
    stop("x must be an object of class 'bn.kcv.list'.")
  # check the connect flag.
  check.logical(connect)

  # store all cross-validation objects in a list for easy handling.
  arg.list = c(list(x), list(...))

  # extract relevant quantities from the cross-validation objects.
  means = numeric(0)
  labels = character(0)
  losses = character(0)

  for (i in seq_along(arg.list)) {

    # transform bn.kcv objects into bn.kcv.list objects.
    if (is(arg.list[[i]], "bn.kcv"))
      arg.list[[i]] = structure(list(arg.list[[i]]), class = "bn.kcv.list")
    else if (!is(arg.list[[i]], "bn.kcv.list"))
      stop("optional argument number", i,
        "must be an object of class 'bn.kcv' or 'bn.kcv.list'.")

    m = sapply(arg.list[[i]], function(x) attr(x, "mean"))
    l = sapply(arg.list[[i]], function(x) attr(x, "loss"))

    means = c(means, m)
    labels = c(labels, rep(as.character(i), length(m)))
    losses = c(losses, l)

  }#FOR

  lattice.cv.bwplot(means = means, labels = labels, losses = losses,
    main = main, xlab = xlab, ylab = ylab, connect = connect)

  invisible(NULL)

}#PLOT.BN.KCV.LIST

