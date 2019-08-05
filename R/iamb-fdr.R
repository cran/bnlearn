
incremental.association.fdr = function(x, cluster = NULL, whitelist,
  blacklist, test, alpha, B, max.sx = ncol(x), complete, debug = FALSE) {

  nodes = names(x)

  # 1. [Compute Markov Blankets]
  mb = smartSapply(cluster, as.list(nodes), ia.fdr.markov.blanket, data = x,
         nodes = nodes, alpha = alpha, B = B, whitelist = whitelist,
         blacklist = blacklist, test = test, max.sx = max.sx,
         complete = complete, debug = debug)
  names(mb) = nodes

  # check markov blankets for consistency.
  mb = bn.recovery(mb, nodes = nodes, mb = TRUE, debug = debug)

  # 2. [Compute Graph Structure]
  mb = smartSapply(cluster, as.list(nodes), neighbour, mb = mb, data = x,
         alpha = alpha, B = B, whitelist = whitelist, blacklist = blacklist,
         test = test, max.sx = max.sx, complete = complete, debug = debug)
  names(mb) = nodes

  # check neighbourhood sets for consistency.
  mb = bn.recovery(mb, nodes = nodes, debug = debug)

  return(mb)

}#INCREMENTAL.ASSOCIATION.FDR

ia.fdr.markov.blanket = function(x, data, nodes, alpha, B, whitelist, blacklist,
  start = character(0), test, max.sx = ncol(x), complete, debug = FALSE) {

  nodes = nodes[nodes != x]
  fdr.threshold = length(nodes) / seq_along(nodes) * sum(1 / seq_along(nodes))
  culprit = character(0)
  whitelisted = nodes[sapply(nodes,
          function(y) { is.whitelisted(whitelist, c(x, y), either = TRUE) })]
  mb = start
  loop.counter = 0
  state = vector(5 * length(nodes), mode = "list")
  last.added = last.removed = NULL

  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* learning the markov blanket of", x, ".\n")

    if (length(start) > 0)
     cat("* initial set includes '", mb, "'.\n")

  }#THEN

  # whitelisted nodes are included by default (if there's a direct arc
  # between them of course they are in each other's markov blanket).
  # arc direction is irrelevant here.
  mb = union(mb, whitelisted)
  # blacklist is not checked, not all nodes in a markov blanket must be
  # neighbours.

  repeat {

    # stop when reaching the maximum size of the conditioning set.
    if (length(mb) > max.sx) {

       if (debug)
         cat("  @ limiting conditioning sets to", max.sx, "nodes.\n")

      break

    }#THEN

    if (ia.detect.infinite.loop(mb, state, loop.counter, debug)) {

      if (!is.null(last.removed)) {

         mb = c(mb, last.removed)
         culprit = c(culprit, last.removed)

      }#THEN
      else if (!is.null(last.added)) {

        mb = setdiff(mb, last.added)
        culprit = c(culprit, last.added)

      }#ELSE

      if (debug)
        cat("  ! ignoring nodes '", culprit, "' from now on.\n")

      # reset the state list so that no further errors are raised.
      state[[loop.counter]] = NULL
      # reset the loop counter to match.
      loop.counter = loop.counter - 1

      warning("prevented infinite loop in Markov blanket learning (node '", x, "').")

    }#THEN

    # increment the loop counter.
    loop.counter = loop.counter + 1

    # save the current markov blanket to detect changes and avoid infinite loops.
    state[[loop.counter]] = mb

    # get an association measure for each of the available nodes.
    association = sapply(nodes, function(node) {
       indep.test(x, node, sx = setdiff(mb, node), data = data, test = test,
                        B = B, alpha = alpha, complete = complete)})
    names(association) = nodes

    # sort the p-values and the FDR thresholds.
    association = association[order(association)]
    names(fdr.threshold) = names(association)

    if (debug) {

      cat("  * computing and sorting p-values.\n")
      sapply(names(association),
        function(x) {
          cat("    >", x, "has p-value", association[x], "with threshold",
            alpha / fdr.threshold[x], ".\n")
        })

    }#THEN

    # remove nodes from the markov blanket (excluding whitelisted nodes and the
    # node added in the last iteration) in order of increasing association.
    if (debug)
      cat("  * checking nodes for exclusion.\n")

    candidates = setdiff(mb, c(whitelisted, last.added, culprit))

    for (node in rev(intersect(names(association), candidates))) {

      if (association[node] * fdr.threshold[node] > alpha) {

        if (debug)
          cat("    @", node, "removed from the markov blanket.\n")

        mb = setdiff(mb, node)
        last.added = NULL
        last.removed = node
        break

      }#THEN
      else {

        if (debug)
          cat("    >", node, "remains in the markov blanket.\n")

      }#ELSE

    }#FOR

    # start again from the top if the markov blanket has changed.
    if (!identical(mb, state[[loop.counter]]))
      next

    # add nodes to the markov blanket in order of decreasing association.
    if (debug)
      cat("  * checking nodes for association.\n")

    candidates = setdiff(nodes, c(mb, last.removed, culprit))

    for (node in intersect(names(association), candidates)) {

      if (association[node] * fdr.threshold[node] <= alpha) {

        if (debug)
         cat("    @", node, "added to the markov blanket.\n")

        mb = c(mb, node)
        last.added = node
        last.removed = NULL
        break

      }#THEN
      else {

        if (debug)
         cat("    >", node, "not added to the markov blanket.\n")

      }#ELSE

    }#FOR

    # if the markov blanket is unchanged, learning is complete.
    if (identical(mb, state[[loop.counter]]))
      break

  }#REPEAT

  return(mb)

}#IA.FDR.MARKOV.BLANKET

ia.detect.infinite.loop = function(mb, state, loop.counter, debug) {

  for (prev.mb in state[seq_len(loop.counter)]) {

    if (!setequal(mb, prev.mb))
      next

    if (debug) {

      cat("  ! recurring markov blanket configuration detected (", mb, ").\n")
      cat("  ! retracing the steps of the learning process:\n")
      sapply(state[seq(loop.counter)],
        function(str) {
          cat("    >", paste(str, collapse = " "), "\n")
        })
     cat("    >", paste(mb, collapse = " "), "\n")

    }#THEN

    return(TRUE)

  }#FOR

  return(FALSE)

}#IA.FDR.DETECT.INFINITE.LOOP

