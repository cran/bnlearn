
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

