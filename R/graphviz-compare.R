
# graphical comparison of different network structures.
graphviz.compare.backend = function(netlist, nodes, groups, layout, shape, main,
    sub, diff, diff.args) {

  # merge and deduplicate all the arcs.
  arclist = lapply(netlist, function(net) {

    if (is(net, "bn"))
      return(net$arcs)
    else
      return(net)

  })

  merged = empty.graph(nodes)
  arcs(merged, check.cycles = FALSE) =
    unique.arcs(do.call("rbind", arclist), nodes = nodes)

  # lay out the graph.
  gr = graphviz.backend(nodes, merged$arcs, groups = groups, layout = layout,
         shape = shape, render = FALSE)
  # extract the labels of the arcs of the first network.
  ref.arcs = apply(arclist[[1]], 1, paste, collapse = "~")
  # extract the labels of the arcs of the merged network.
  grlabels = names(graph::edgeRenderInfo(gr)[["splines"]])
  graph::edgeRenderInfo(gr)[["lty"]] = "solid"

  # allocate the return value.
  graphlist = vector(length(arclist), mode = "list")

  # iterate over the network structures.
  for (i in seq_along(arclist)) {

    # keep the original and modify a copy.
    gr.temp = gr
    edges.temp = graph::edgeRenderInfo(gr.temp)
    # extract the labels of the arcs of the current network.
    cur.arcs = arcs2grlabels(arclist[[i]])

    # make sure arrowheads appear in the right places; merging two directed
    # arcs, or one directed and one undirected arc, can make the direction
    # disappear.
    directed = which.directed(arclist[[i]], nodes)
    und.arcs = intersect(cur.arcs[!directed], grlabels)

    edges.temp[["arrowhead"]][und.arcs] = "none"
    edges.temp[["arrowtail"]][und.arcs] = "none"

    dir.arcs = grlabels[grlabels %in%
                 arcs2grlabels(arclist[[i]][directed, , drop = FALSE], both = TRUE)]

    for (arc in dir.arcs) {

      if (arc %in% cur.arcs) {

        edges.temp[["direction"]][arc] = "forward"
        edges.temp[["arrowhead"]][arc] = "open"
        edges.temp[["arrowtail"]][arc] = "none"

      }#THEN
      else {

        edges.temp[["direction"]][arc] = "back"
        edges.temp[["arrowhead"]][arc] = "none"
        edges.temp[["arrowtail"]][arc] = "open"

        # this complicated dance with the order of the spline segments is
        # needed to work around a bug in Rgraphviz, which implictly assumes
        # that all the points are in a single spline segment and thus botches
        # the placement of the arrow when the direction is "back".
        splines = edges.temp[["splines"]][[arc]]
        splines = rev(splines)
        edges.temp[["splines"]][[arc]] = splines

      }#ELSE

    }#FOR

    # go ahead with the formatting.
    if (diff == "none") {

      # no formatting, just plot all the network structures.
      edges.temp[["col"]][] = "transparent"
      edges.temp[["col"]][und.arcs] = "black"
      edges.temp[["col"]][dir.arcs] = "black"

    }#THEN
    else if (diff == "from-first") {

      if (i == 1) {

        edges.temp[["col"]][] = "transparent"
        edges.temp[["col"]][und.arcs] = "black"
        edges.temp[["col"]][dir.arcs] = "black"

      }#THEN
      else {

        edges.temp[["col"]][] = "transparent"

        # classify arcs as true positives, true negatives and false negatives.
        sorted = compare.backend(arclist[[1]], arclist[[i]], nodes, arcs = TRUE)

        # remove false negatives that have a matching false positive, formatting
        # is easier if the two sets are disjoint.
        dupes = sorted$fp[which.listed(sorted$fp, sorted$fn, either = TRUE), , drop = FALSE]
        sorted$fn = sorted$fn[!which.listed(sorted$fn, dupes, either = TRUE), , drop = FALSE]

        # matching true positive directed arcs.
        tp.labels = grlabels[grlabels %in% arcs2grlabels(sorted$tp, both = TRUE)]

        if (length(tp.labels) > 0) {

          edges.temp[["col"]][tp.labels] = diff.args$tp.col
          edges.temp[["lty"]][tp.labels] = diff.args$tp.lty
          if ("tp.lwd" %in% names(diff.args))
            edges.temp[["lwd"]][tp.labels] = diff.args$tp.lwd

        }#THEN

        # match false positive directed arcs.
        fp.labels = grlabels[grlabels %in% arcs2grlabels(sorted$fp, both = TRUE)]

        if (length(fp.labels) > 0) {

          edges.temp[["col"]][fp.labels] = diff.args$fp.col
          edges.temp[["lty"]][fp.labels] = diff.args$fp.lty
          if ("fp.lwd" %in% names(diff.args))
            edges.temp[["lwd"]][fp.labels] = diff.args$fp.lwd

        }#THEN

        # match false negative directed arcs.
        fn.labels.fwd = grlabels[grlabels %in% arcs2grlabels(sorted$fn)]
        fn.labels.bwd = grlabels[grlabels %in% arcs2grlabels(sorted$fn[, 2:1, drop = FALSE])]

        if (length(c(fn.labels.fwd, fn.labels.bwd)) > 0) {

          edges.temp[["col"]][c(fn.labels.fwd, fn.labels.bwd)] = diff.args$fn.col
          edges.temp[["lty"]][c(fn.labels.fwd, fn.labels.bwd)] = diff.args$fn.lty
          if ("fn.lwd" %in% names(diff.args))
            edges.temp[["lwd"]][c(fn.labels.fwd, fn.labels.bwd)] = diff.args$fn.lwd

        }#THEN

        # arcs that do not appear in the current graph should be directed as in
        # the reference graph, if they are not.
        both = names(which(edges.temp[["direction"]] == "both"))

        for (arc in both) {

          if (arc %in% setdiff(fn.labels.fwd, fn.labels.bwd)) {

            edges.temp[["direction"]][arc] = "forward"
            edges.temp[["arrowhead"]][arc] = "open"
            edges.temp[["arrowtail"]][arc] = "none"

          }#THEN
          else if (arc %in% setdiff(fn.labels.bwd, fn.labels.fwd)) {

            edges.temp[["direction"]][arc] = "back"
            edges.temp[["arrowhead"]][arc] = "none"
            edges.temp[["arrowtail"]][arc] = "open"
            # same complicated dance as above.
            splines = edges.temp[["splines"]][[arc]]
            splines = rev(splines)
            edges.temp[["splines"]][[arc]] = splines

          }#THEN

        }#FOR

      }#ELSE

    }#THEN

    graph::edgeRenderInfo(gr.temp) = edges.temp

    # add the title and subtitle, to label the panels.
    graph::graphRenderInfo(gr.temp)$main = main[i]
    graph::graphRenderInfo(gr.temp)$sub = sub[i]

    # save the formatted network to return it later.
    graphlist[[i]] = gr.temp
    # plot the formatted network structure, unless told not to do that.
    if ((i == 1) && (diff == "from-first") && !diff.args$show.first)
      next

    Rgraphviz::renderGraph(gr.temp)

  }#FOR

  invisible(graphlist)

}#GRAPHVIZ.COMPARE.BACKEND

# utility function to paster graphviz labels together.
arcs2grlabels = function(arcset, both = FALSE) {

  if (both)
    apply(arcs.rbind(arcset, arcset, reverse2 = TRUE), 1, paste, collapse = "~")
  else
    apply(arcset, 1, paste, collapse = "~")

}#ARCS2GRLABEL
