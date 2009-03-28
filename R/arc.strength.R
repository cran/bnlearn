
arc.strength.test = function(network, data, test, alpha, debug) {

  drop = function(arc) {

    if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("* computing strength for arc", arc[1], "->", arc[2], ".\n")

    }#THEN

    parents =
      network$nodes[[arc[2]]]$parents[network$nodes[[arc[2]]]$parents != arc[1]]

    a = conditional.test(arc[1], arc[2], parents, data = data, test = test)

    if (debug) {

      cat("  > testing", arc[1], "->", arc[2],
        "with conditioning set '", parents, "'.\n")
      cat("    > p-value is", a, ".\n")

    }#THEN

    return(a)

  }#DROP

  if (debug) {

    cat("----------------------------------------------------------------\n")
    print(network)

  }#THEN

  # populate the strength data frame.
  strength = data.frame(network$arcs, strength = apply(network$arcs, 1, drop),
               stringsAsFactors = FALSE)

  return(strength)

}#ARC.STRENGTH.TEST

arc.strength.score = function(network, data, score, extra, debug) {

  drop = function(arc) {

    if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("* computing strength for arc", arc[1], "->", arc[2], ".\n")

    }#THEN

    better = score.delta(arc = arc, network = network, data = data,
               score = score, score.delta = 0,
               reference.score = reference.score, op = "drop",
               extra = extra)

    if (debug) {

      cat("  > updated score for node", arc[2], "is", better$updates, ".\n")
      cat("  > score delta", better$delta, ".\n")

    }#THEN

    return(better$delta)

  }#DROP

  # cache nodes' labels.
  nodes = names(data)
  # set the reference score.
  reference.score = per.node.score(network = network, score = score,
                      nodes = nodes, extra.args = extra, data = data)

  if (debug) {

    cat("----------------------------------------------------------------\n")
    print(network)
    cat("* current score:", sum(reference.score), "\n")
    cat("----------------------------------------------------------------\n")
    cat("* original scores of the nodes of the graphs are:\n")
    for (n in nodes)
      cat("  > original score for node", n, "is", reference.score[n], ".\n")

  }#THEN

  # populate the strength data frame.
  strength = data.frame(network$arcs, strength = apply(network$arcs, 1, drop),
               stringsAsFactors = FALSE)

  return(strength)

}#ARC.STRENGTH.SCORE

# convert an arc strength object to the corresponding line widths for plotting.
strength2lwd = function(strength, threshold, cutpoints, debug = TRUE) {

  s = strength[, "strength"]
  mode = attr(strength, "mode")

  # sanitize user-defined cut points, if any.
  if (!missing(cutpoints)) {

    if (!is.numeric(cutpoints) || any(is.nan(cutpoints)))
      stop("cut points must be numerical values.")
    if (length(s) <= length(cutpoints))
      stop("there are at least as many cut points as strength values.")

  }#THEN

  if (debug) {

    cat("* using threshold:", threshold, "\n")
    cat("* reported arc strength are:\n")
    print(strength)

  }#THEN

  if (mode == "test") {

    # use user-defined cut points if available.
    if (missing(cutpoints))
      cutpoints = unique(c(0, threshold/c(10, 5, 2, 1.5, 1), 1))
    else
      cutpoints = sort(cutpoints)

    # p-values are already scaled, so the raw quantiles are good cut points.
    arc.weights = cut(s, cutpoints, labels = FALSE, include.lowest = TRUE)

    arc.weights = length(cutpoints) - arc.weights

  }#THEN
  else if (mode == "score") {

    # score deltas are defined on a reversed scale (the more negative
    # the better); change their sign (and that of the threshold)
    # for simplicity
    threshold = -threshold
    s = -s

    # define a set of cut points using the quantiles from the empirical
    # distribution of the negative score deltas (that is, the ones
    # corresponding to significant arcs) or use user-defined ones.
    if (missing(cutpoints)) {

      significant = s[s > threshold]
      q = quantile(significant, c(0.50, 0.75, 0.90, 0.95, 1), names = FALSE)
      cutpoints = sort(c(-Inf, threshold, unique(q), Inf))

    }#THEN

    arc.weights = cut(s, cutpoints, labels = FALSE)

  }#THEN

  # arcs beyond the significance threshold are given a negative weight,
  # so that graphviz.backend() will draw them as dashed lines.
  arc.weights[arc.weights == 1] = -1

  if (debug) {

    cat("* using cut points for strength intervals:\n")
    print(cutpoints)
    cat("* arc weights:", arc.weights, "\n")

  }#THEN

  return(arc.weights)

}#STRENGTH2LWD
