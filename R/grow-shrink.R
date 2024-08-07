
grow.shrink = function(x, cluster = NULL, whitelist, blacklist, test, alpha,
    extra.args = list(), max.sx = ncol(x), complete, debug = FALSE) {

  nodes = names(x)

  # 1. [Compute Markov Blankets]
  mb = smartSapply(cluster, as.list(nodes), gs.markov.blanket, data = x,
         nodes = nodes, alpha = alpha, extra.args = extra.args,
         whitelist = whitelist, blacklist = blacklist, test = test,
         max.sx = max.sx, debug = debug)
  names(mb) = nodes

  # check markov blankets for consistency.
  mb = bn.recovery(mb, nodes = nodes, mb = TRUE, debug = debug)

  # 2. [Compute Graph Structure]
  mb = smartSapply(cluster, as.list(nodes), neighbour, mb = mb, data = x,
         alpha = alpha, extra.args = extra.args, whitelist = whitelist,
         blacklist = blacklist, test = test, max.sx = max.sx, debug = debug)
  names(mb) = nodes

  # check neighbourhood sets for consistency.
  mb = bn.recovery(mb, nodes = nodes, debug = debug)

  return(mb)

}#GROW.SHRINK

gs.markov.blanket = function(x, data, nodes, alpha, extra.args = list(),
    whitelist, blacklist, start = character(0), test, max.sx = ncol(x),
    debug = FALSE) {

  nodes = nodes[nodes != x]
  whitelisted = nodes[sapply(nodes,
          function(y) { is.whitelisted(whitelist, c(x, y), either = TRUE) })]
  mb = start

  # growing phase
  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* learning the markov blanket of", x, ".\n")

    if (length(start) > 0)
      cat("* initial set includes '", mb, "'.\n")

  }#THEN

  # whitelisted nodes are included by default (if there's a direct arc
  # between them of course they are in each other's markov blanket).
  # arc direction is irrelevant here.
  mb = unique(c(mb, whitelisted))
  nodes = nodes[nodes %!in% mb]
  # blacklist is not checked, not all nodes in a markov blanket must be
  # neighbours.

  # grow phase.
  repeat {

    # store the current size of the Markov Blanket.
    mb.size = length(mb)

    for (y in nodes) {

      # avoid testing with large conditioning sets.
      if (mb.size > max.sx) {

        if (debug)
          cat("  * skipping node", y,
              "because the tests would involve a conditioning set larger than",
              max.sx, ".\n")

        next

      }#THEN

      if (debug)
        cat("  * checking node", y, "for inclusion.\n")


      a = indep.test(x, y, mb, data = data, test = test,
            extra.args = extra.args, alpha = alpha)

      if (a <= alpha) {

        # add the node to the Markov blanket.
        mb = c(mb, y)
        # do not check the same node again.
        nodes = nodes[nodes != y]

        if (debug) {

          cat("    > node", y, "included in the markov blanket ( p-value:", a, ").\n")
          cat("    > markov blanket (", length(mb), "nodes ) now is '", mb, "'.\n")
          cat("    > restarting grow loop.\n")

        }#THEN

        break

      }#THEN
      else if (debug) {

        cat("    >", x, "indep.", y, "given '", mb, "' ( p-value:", a, ").\n")

      }#THEN

    }#FOR

    # if the Markov blanket is unchanged exit the grow phase.
    if (length(mb) == mb.size) break

  }#REPEAT

  # whitelisted nodes are neighbours, they cannot be removed from the
  # markov blanket.
  nodes = mb[mb %!in% whitelisted]

  # shrinking phase.
  repeat {

    # store the current size of the Markov Blanket.
    mb.size = length(mb)

    for (y in nodes) {

      if (debug)
        cat("  * checking node", y, "for exclusion (shrinking phase).\n")

      a = indep.test(x, y, mb[mb != y], data = data, test = test,
            extra.args = extra.args, alpha = alpha)

      if (a > alpha) {

        # update the markov blanket.
        mb = mb[mb != y]
        # do not check the same node again.
        nodes = nodes[nodes != y]

        if (debug) {

          cat("    > node", y, "removed from the markov blanket. ( p-value:", a, ")\n")
          cat("    > conditioning subset: '", mb, "'\n")
          cat("    > restarting shrink loop.\n")

        }#THEN

        break

      }#THEN
      else if (debug) {

        cat("    > node", y, "remains in the markov blanket. ( p-value:", a, ")\n")

      }#THEN

    }#FOR

    # if the Markov blanket is unchanged exit the grow phase.
    if (length(mb) == mb.size) break

  }#REPEAT

  return(mb)

}#GS.MARKOV.BLANKET

