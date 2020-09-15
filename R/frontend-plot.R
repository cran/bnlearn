
# generic plot of an object of class 'bn' or 'bn.fit' using graphviz.
graphviz.plot = function(x, highlight = NULL, groups, layout = "dot",
    shape = "circle", main = NULL, sub = NULL, render = TRUE) {

  # check whether graphviz is loaded.
  check.and.load.package("Rgraphviz")
  # check x's class.
  check.bn.or.fit(x)
  # check render.
  check.logical(render)

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
    groups = groups, layout = layout, shape = shape, main = main, sub = sub,
    render = render)

}#GRAPHVIZ.PLOT

# plot a graph with arcs formatted according to their own strength.
strength.plot = function(x, strength, threshold, cutpoints, highlight = NULL,
    groups, layout = "dot", shape = "circle", main = NULL, sub = NULL,
    render = TRUE, debug = FALSE) {

  # check whether graphviz is loaded.
  check.and.load.package("Rgraphviz")
  # check x's class.
  check.bn(x)
  # check the strength parameter.
  check.bn.strength(strength, nodes = names(x$nodes))
  # check the strength threshold.
  threshold = check.threshold(threshold, strength)

  # check and match the strength coefficients.
  str = match.arcs.and.strengths(arcs = x$arcs, nodes = names(x$nodes),
          strengths = strength)

  # sanitize user-defined cut points, if any.
  method =  attr(strength, "method")
  cutpoints = check.cutpoints(cutpoints, strength = str, threshold = threshold,
                method = method)

  # compute arc weights from the bn.strength object.
  arc.weights = strength2lwd(strength = str, threshold = threshold,
                  cutpoints = cutpoints, method = method,
                  arcs = x$arcs, debug = debug)

  # create and maybe plot the graph object.
  gr = graphviz.backend(nodes = names(x$nodes), arcs = x$arcs,
         highlight = highlight, groups = groups, arc.weights = arc.weights,
         layout = layout, shape = shape, main = main, sub = sub,
         render = render)

  # save the arc strengths in the weights.
  graph::edgeData(gr, from = x$arcs[, 1], to = x$arcs[, 2],
                     attr = "weight") = str

  # return the graph object for further customization.
  invisible(gr)

}#STRENGTH.PLOT

# combine barcharts and graph displays in a single plot.
graphviz.chart = function(x, type = "barchart", layout = "dot",
    draw.levels = TRUE, grid = FALSE, scale = c(0.75, 1.1), col = "black",
    bg = "transparent", text.col = "black", bar.col = "black", strip.bg = bg,
    main = NULL, sub = NULL) {

  nodes = nodes(x)
  nnodes = length(nodes)

  # check whether graphviz is loaded.
  check.and.load.package("Rgraphviz")
  # check whether gRain is loaded.
  check.and.load.package("gRain")
  # check x's class.
  check.fit(x)
  # only discrete networks are supported.
  if (!is(x, c("bn.fit.dnet", "bn.fit.onet", "bn.fit.donet")))
    stop("only discrete networks are supported.")

  # check the plot type.
  check.label(type, choices = c("barchart", "dotplot", "barprob"),
    argname = "node display type", see = "graphviz.chart")
  # check the levels switch.
  check.logical(draw.levels)
  # check the grid.
  grid = check.quantile.grid(grid)
  # check the nodes scale factors.
  if (!is.positive.vector(scale) || (length(scale) != 2))
    stop("the horizontal and vertical scale factors for the nodes must be two positive numbers.")

  # check the colours.
  col = check.colour(col, num = nnodes, expand = TRUE, labels = nodes)
  bg = check.colour(bg, num = nnodes, expand = TRUE, labels = nodes)
  text.col = check.colour(text.col, num = nnodes, expand = TRUE, labels = nodes)
  bar.col = check.colour(bar.col, num = nnodes, expand = TRUE, labels = nodes)
  strip.bg = check.colour(strip.bg, num = nnodes, expand = TRUE, labels = nodes)

  graphviz.chart.backend(fitted = x, type = type, layout = layout,
    draw.levels = draw.levels, grid = grid, scale = scale, col = col, bg = bg,
    text.col = text.col, bar.col = bar.col, strip.bg = strip.bg, main = main,
    sub = sub)

}#GRAPHVIZ.CHART

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

plot.bn.fit = function(x, ...) {

  stop("no plot() method available for 'bn.fit' objects, see ?`bn.fit plots` for alternatives.")

}#PLOT.BN.FIT

# ECDF plot for 'bn.strength' objects.
plot.bn.strength = function(x, draw.threshold = TRUE, main = NULL,
    xlab = "arc strengths", ylab = "CDF(arc strengths)", ...) {

  # check the draw.threshold flag.
  check.logical(draw.threshold)
  # check the arcs strengths.
  check.bn.strength(x, valid = c("bootstrap", "bayes-factor"))

  # match the dots arguments.
  dots = sanitize.plot.dots(eval(list(...)), meaningless = c("xlim", "ylim"))

  # handle the title.
  if (is.null(main))
    main = paste("threshold = ", round(attr(x, "threshold"), digits = 3))

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

# graphical comparison of different network structures.
graphviz.compare = function(x, ..., groups, layout = "dot", shape = "circle",
    main = NULL, sub = NULL, diff = "from-first", diff.args = list()) {

  available.diff.methods = c("none", "from-first")

  # check whether graphviz is loaded.
  check.and.load.package("Rgraphviz")
  # check the first network structure.
  check.bn(x)
  nodes = names(x$nodes)
  # collect all the networks in a singe list.
  netlist = c(list(x), list(...))
  # check that the networks in the list are valid and agree with each other.
  check.customlist(netlist, nodes = nodes)
  # check the titles and the subtitles for the networks.
  if (!is.null(main))
    if (!is.string.vector(main) || (length(main) != length(netlist)))
      stop("'main' must a vector of character strings, one for each network.")
  if (!is.null(sub))
    if (!is.string.vector(sub) || (length(sub) != length(netlist)))
      stop("'sub' must a vector of character strings, one for each network.")
  # check the diff method.
  check.label(diff, choices = available.diff.methods, argname = "diff")
  # check the diff extra arguments.
  if (diff == "none") {

    # nothing to do.
    check.unused.args(diff.args, character(0))

  }#THEN
  else if (diff == "from-first") {

    args = c("tp.col", "tp.lty", "tp.lwd", "fp.col", "fp.lty", "fp.lwd",
             "fn.col", "fn.lty", "fn.lwd", "show.first")

    # check unused extra arguments.
    check.unused.args(diff.args, args)

    # check the line type of the different classes of arcs.
    if ("tp.lty" %in% names(diff.args))
      check.lty(diff.args$tp.lty)
    else
      diff.args$tp.lty = "solid"
    if ("fp.lty" %in% names(diff.args))
      check.lty(diff.args$fp.lty)
    else
      diff.args$fp.lty = "solid"
    if ("fn.lty" %in% names(diff.args))
      check.lty(diff.args$fn.lty)
    else
      diff.args$fn.lty = "dashed"

    # check the colour of the different classes of arcs.
    if ("tp.col" %in% names(diff.args))
      check.colour(diff.args$tp.col)
    else
      diff.args$tp.col = "black"

    if ("fp.col" %in% names(diff.args))
      check.colour(diff.args$fp.col)
    else
      diff.args$fp.col = "red"

    if ("fn.col" %in% names(diff.args))
      check.colour(diff.args$fn.col)
    else
      diff.args$fn.col = "blue"

    # check the line width of the different classes of arcs.
    for (lwd in c("tp.lwd", "fp.lwd", "fn.lwd"))
      if (lwd %in% names(diff.args))
        if (!is.positive(diff.args[[lwd]]))
          stop("diff.args$", lwd, " must be a positive number.")

    # check whether the first, reference network should be plotted.
    if ("show.first" %in% names(diff.args))
      check.logical(diff.args$show.first)
    else
      diff.args$show.first = TRUE

  }#THEN

  # the sanitization of "layout" and "shape" is left to the backend.

  graphviz.compare.backend(netlist = netlist, nodes = nodes, groups = groups,
    layout = layout, shape = shape, main = main, sub = sub, diff = diff,
    diff.args = diff.args)

}#GRAPHVIZ.COMPARE

