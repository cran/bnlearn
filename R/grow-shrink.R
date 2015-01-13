
grow.shrink.optimized = function(x, whitelist, blacklist, test, alpha,
  B, strict, debug = FALSE) {

  nodes = names(x)
  mb2 = mb = list()

  # 1. [Compute Markov Blankets]
  for (node in nodes) {

    backtracking = unlist(sapply(mb, function(x){ node %in% x  }))

    mb[[node]] = gs.markov.blanket(node, data = x, nodes = nodes,
         alpha = alpha, B = B, whitelist = whitelist, blacklist = blacklist,
         backtracking = backtracking, test = test, debug = debug)

  }#FOR

  # check markov blankets for consistency.
  mb = bn.recovery(mb, nodes = nodes, strict = strict, mb = TRUE, debug = debug)

  # 2. [Compute Graph Structure]
  for (node in nodes) {

    backtracking = unlist(sapply(mb2, function(x){ node %in% x[["nbr"]]  }))

    # save results in a copy of mb.
    mb2[[node]] = neighbour(node, mb = mb, data = x, alpha = alpha,
         B = B, whitelist = whitelist, blacklist = blacklist,
         backtracking = backtracking, test = test, debug = debug)

  }#FOR

  # update mb with the results of neighbour().
  mb = mb2

  # check neighbourhood sets for consistency.
  mb = bn.recovery(mb, nodes = nodes, strict = strict, debug = debug)

  return(mb)

}#GROW.SHRINK.OPTIMIZED

grow.shrink = function(x, cluster = NULL, whitelist, blacklist, test, alpha, B,
  strict, debug = FALSE) {

  nodes = names(x)

  # 1. [Compute Markov Blankets]
  mb = smartLapply(cluster, as.list(nodes), gs.markov.blanket, data = x,
         nodes = nodes, alpha = alpha, B = B, whitelist = whitelist,
         blacklist = blacklist, test = test, debug = debug)
  names(mb) = nodes

  # check markov blankets for consistency.
  mb = bn.recovery(mb, nodes = nodes, strict = strict, mb = TRUE, debug = debug)

  # 2. [Compute Graph Structure]
  mb = smartLapply(cluster, as.list(nodes), neighbour, mb = mb, data = x,
         alpha = alpha, B = B, whitelist = whitelist, blacklist = blacklist,
         test = test, debug = debug)
  names(mb) = nodes

  # check neighbourhood sets for consistency.
  mb = bn.recovery(mb, nodes = nodes, strict = strict, debug = debug)

  return(mb)

}#GROW.SHRINK

gs.markov.blanket = function(x, data, nodes, alpha, B, whitelist, blacklist,
  start = character(0), backtracking = NULL, test, debug = FALSE) {

  nodes = nodes[nodes != x]
  known.good = known.bad = c()
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

  # use backtracking for a further screening of the nodes to be checked.
  if (!is.null(backtracking)) {

    # nodes whose markov blanket includes this node are included, because
    # X \in MB(Y) <=> Y \in MB(X); they can be removed in the backward phase
    # later on if they are false positives.
    known.good = names(backtracking[backtracking])
    mb = unique(c(mb, known.good))

    # known.good nodes are not to be checked for inclusion, and "mb" is used
    # for the shrinking phase instead of "nodes", we can just remove them.
    nodes = nodes[nodes %!in% known.good]

    # and vice versa X \not\in MB(Y) <=> Y \not\in MB(X)
    known.bad = names(backtracking[!backtracking])

    if (debug) {

      cat("    * known good (backtracking): '", known.good, "'.\n")
      cat("    * known bad (backtracking): '", known.bad, "'.\n")
      cat("    * nodes still to be tested for inclusion: '", nodes, "'.\n")

    }#THEN

  }#THEN

  # grow phase.
  repeat {

    # store the current size of the Markov Blanket.
    mb.size = length(mb)

    for (y in nodes) {

      if (debug)
        cat("  * checking node", y, "for inclusion.\n")

      a = indep.test(x, y, mb, data = data, test = test, B = B,
            alpha = alpha)

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
  # markov blanket; on the other hand, known.good nodes from backtracking
  # should be checked to remove false positives.
  nodes = mb[mb %!in% whitelisted]

  # shrinking phase.
  repeat {

    # store the current size of the Markov Blanket.
    mb.size = length(mb)

    for (y in nodes) {

      if (debug)
        cat("  * checking node", y, "for exclusion (shrinking phase).\n")

      a = indep.test(x, y, mb[mb != y], data = data, test = test, B = B,
            alpha = alpha)

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

