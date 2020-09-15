
graphviz.chart.backend = function(fitted, type = "barchart", layout = "dot",
    draw.levels = TRUE, grid = NULL, scale = c(0.75, 1.1), col = "black",
    bg = "transparent", text.col = "black", bar.col = "black",
    strip.bg = "transparent", main = NULL, sub = NULL) {

  # sanitize the graph layout.
  check.label(layout, choices = graphviz.layouts, argname = "graph layout")

  # create the graphNEL object.
  nodes = names(fitted)
  nnodes = length(nodes)
  arcs = fit2arcs(fitted)

  graph.obj = new("graphNEL", nodes = nodes, edgeL = arcs2elist(arcs, nodes),
                edgemode = 'directed')

  # create the marginal probabilities.
  probs = gRain::querygrain(as.grain(fitted), nodes = nodes, type = "marginal")

  # create the function that will draw the charts (the arguments are mandated
  # by Rgraphviz, and are not really used apart from the first one).
  drawFuns = function(node, ur, attrs, radConv) {

      nc = Rgraphviz::getNodeCenter(node)
      nl = node@txtLabel@labelText
      chartGlyph(probs[[nl]], xpos = Rgraphviz::getX(nc),
          ypos = Rgraphviz::getY(nc), node = nl, height = node@height,
          width = node@rWidth + node@lWidth, draw.labels = draw.levels,
          grid = grid, type = type, max.levels = max(sapply(probs, length)),
          col = col[nl], bg = bg[nl], text.col = text.col,
          bar.col = bar.col[nl], strip.bg = strip.bg[nl])

    }#FUNCTION

  # initialize the plot, and compute the margins for the title and the subtitle.
  if (names(dev.cur()) == "null device")
    plot.new()

  mai.sub = mai.title = 0

  if (!is.null(main))
    mai.title = strheight(main, "inches") + 0.2
  if (!is.null(sub))
    mai.sub = strheight(sub, "inches") + 0.2

  mai = c(mai.sub, 0, mai.title, 0)

  # draw the plot.
  Rgraphviz::plot(graph.obj, layout, drawNode = drawFuns, mai = mai,
    nodeAttrs = list(
      shape = structure(rep("rectangle", nnodes), names = nodes),
      height = structure(rep(scale[1], nnodes), names = nodes),
      width = structure(rep(scale[2], nnodes), names = nodes)
    ))

  # add the title and the subtitle; Rgraphviz messes up the placement of the
  # subtitle, so it must be handled separately.
  if (!is.null(main))
    title(main = main, line = 0.5)
  if (!is.null(sub))
    title(sub = sub, line = 0.25)

  invisible(NULL)

}#GRAPHVIZ.CHART.BACKEND

# draw the chart for a single node in a graphviz plot.
chartGlyph = function(prob, xpos, ypos, height, width, node, draw.labels, grid,
    type, max.levels, col, bg, text.col, bar.col, strip.bg) {

  # compute the boundaries of the box and the height of the title box.
  xlim = xpos + c(-1, 1) * width / 2
  ylim = ypos + c(-1, 1) * height / 2
  title.box.height = height * 0.22

  # draw the background, if any.
  if (bg != "transparent")
    rect(xleft = xlim[1], ybottom = ylim[1], xright = xlim[2],
       ytop = ylim[2] - title.box.height, col = bg, border = "transparent")
  if (strip.bg != "transparent")
    rect(xleft = xlim[1], ybottom = ylim[2] - title.box.height,
       xright = xlim[2], ytop = ylim[2], col = strip.bg,
       border = "transparent")

  # place the label of the node at the top of the box, finding the best cex.
  best.cex = largest.cex(node, height = title.box.height, width = width)
  text(x = xpos, y = ylim[2] - title.box.height / 2, node, cex = best.cex,
    col = text.col)

  # move below the box title.
  y.top = ylim[2] - title.box.height
  # compute the vertical distance between the bars/lines.
  stepping = (y.top - ylim[1]) / length(prob)
  # compute the maximum bar height.
  total.bar.height = (y.top - ylim[1]) / max.levels

  if (draw.labels) {

    # extract the labels and set their positions and size.
    labels = names(prob)
    label.box.width = width * 0.30
    label.x = xlim[1] + label.box.width
    # delimit the area in which to draw the bars/lines, making sure it does not
    # overlap either the labels or the bounding box.
    x.left = xlim[1] + label.box.width + 0.03 * width
    total.bar.width = xlim[2] - x.left - 0.03 * width
    # set the font size of the labels (at 80% of bar height in barchart).
    best.cex = min(sapply(labels, largest.cex,
                     height = total.bar.height, hfrac = 0.7 * 0.8,
                     width = label.box.width), wfrac = 0.95)

  }#THEN
  else {

    # delimit the area in which to draw the bars/lines, making sure it does not
    # overlap either the labels or the bounding box.
    x.left = xlim[1] + 0.03 * width
    total.bar.width = xlim[2] - x.left - 0.03 * width

  }#ELSE

  # draw a grid to make probabilities easier to read.
  if (!is.null(grid)) {

    grid.points = qunif(grid, min = x.left, max = x.left + total.bar.width)

    for (g in grid.points)
      lines(x = rep(g, 2), y = c(ylim[1], y.top), col = lighter.colour(col, 0.75))

  }#THEN

  # the left x-coord is fixed, compute the right x-coord and the y-coord.
  x.right = x.left + prob * total.bar.width
  y.pos = y.top - (seq_along(prob)  - 1/2) * stepping

  # draw the bars/lines representing the probabilities.
  if (type == "dotplot") {

    # draw the lines that fill the role of the bars (faster in a for loop).
    for (i in seq_along(prob))
      lines(x = c(x.left, x.right[i]), y = rep(y.pos[i], 2),
        col = bar.col, lwd = 2)
    # add a bullets at the end of each line.
    symbols(x = x.right, y = y.pos, inches = FALSE, bg = bar.col, fg = bar.col,
      circles = rep(total.bar.height * 0.15, length(prob)), add = TRUE)

  }#THEN
  else if (type == "barchart") {

    # draw a rectangle for the bar.
    bar.bottom = y.pos - 0.35 * total.bar.height
    bar.top = y.pos + 0.35 * total.bar.height
    rect(xleft = x.left, ybottom = bar.bottom, xright = x.right, ytop = bar.top,
      col = lighter.colour(bar.col), border = bar.col)

  }#THEN
  else if (type == "barprob") {

    # draw a rectangle for the bar.
    bar.bottom = y.pos - 0.45 * total.bar.height
    bar.top = y.pos + 0.45 * total.bar.height
    rect(xleft = x.left, ybottom = bar.bottom, xright = x.right, ytop = bar.top,
      col = lighter.colour(bar.col, 0.75), border = lighter.colour(bar.col))
    # print the probabilities on top of the bars.
    prob.string = sprintf("%.3f", lrm.round(prob, digits = 3))
    text(x = x.left + 0.5 * total.bar.width, y = y.pos, prob.string,
      cex = largest.cex(prob.string[1], total.bar.height, total.bar.width, hfrac = 0.56),
      col = text.col)

  }#THEN

  # place the label to the left of the bar/line, without using "pos" to make
  # it align properly with the bars/lines and the probabilities in "barprob".
  if (draw.labels)
    text(x = label.x - 0.5 * strwidth(labels, cex = best.cex),
      y = y.pos, labels, cex = best.cex, col = text.col)

  # draw the bounding box.
  rect(xleft = xlim[1], ybottom = ylim[1], xright = xlim[2],
     ytop = ylim[2], col = "transparent", border = col)
  lines(x = c(xlim[1], xlim[2]), y = rep(ylim[2] - title.box.height, 2),
     col = col)

  invisible(NULL)

}#CHARTGLYPH

