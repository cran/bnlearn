
inter.incremental.association = function(x, cluster = NULL, whitelist,
  blacklist, test, alpha, B, max.sx = ncol(x), complete, debug = FALSE) {

  nodes = names(x)

  # 1. [Compute Markov Blankets]
  mb = smartSapply(cluster, as.list(nodes), inter.ia.markov.blanket, data = x,
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

}#INTER.INCREMENTAL.ASSOCIATION

inter.ia.markov.blanket = function(x, data, nodes, alpha, B, whitelist,
  blacklist, start = character(0), test, max.sx = ncol(x), complete,
  debug = FALSE) {

  nodes = nodes[nodes != x]
  culprit = character(0)
  whitelisted = nodes[sapply(nodes,
          function(y) { is.whitelisted(whitelist, c(x, y), either = TRUE) })]
  mb = start
  loop.counter = 1
  state = vector(5 * length(nodes), mode = "list")

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

  repeat {

    # stop if there are no nodes left.
    if (length(nodes[nodes %!in% c(mb, culprit)]) == 0 || is.null(nodes))
      break

    if (length(mb) > max.sx) {

       if (debug)
         cat("  @ limiting conditioning sets to", max.sx, "nodes.\n")

      break

    }#THEN

    # get a snapshot of the markov blanket status.
    mb.snapshot = mb

    # get an association measure for each of the available nodes.
    association = indep.test(nodes[nodes %!in% c(mb, culprit)], x, sx = mb,
                    test = test, data = data, B = B, alpha = alpha,
                    complete = complete)

    # stop if there are no candidates for inclusion; the markov blanket
    # would obviously be unchanged.
    if (all(association > alpha))
      break
    # get the one which maximizes the association measure.
    to.add = names(which.min(association))

    if (debug) {

      cat("  * checking nodes for association.\n")
      sapply(names(association),
        function(x) {  cat("    >", x, "has p-value", association[x], ".\n")})
      cat("    @", to.add, "included in the markov blanket ( p-value:",
        association[to.add], ").\n")
      cat("    > markov blanket (", length(mb) + 1, " nodes ) now is '", c(mb, to.add), "'.\n")

    }#THEN

    if (association[to.add] <= alpha) mb = c(mb, to.add)

    if (debug)
      cat("  * checking nodes for exclusion.\n")

    # whitelisted nodes are neighbours, they cannot be removed from the markov
    # blanket; the last node added in phase I will never be removed, because
    # the tests for inclusion and removal are identical.
    fixed = c(to.add, whitelisted)
    fixed = fixed[fixed != ""]

    pv = roundrobin.test(x = x, z = mb, fixed = fixed, data = data, test = test,
           B = B, alpha = alpha, complete = complete, debug = debug)

    mb = intersect(mb, c(names(pv[pv < alpha]), fixed))

    # if the markov blanket did not change, stop iterating.
    if (identical(mb, mb.snapshot)) {

      if (debug)
        cat("  @ markov blanket is unchanged, stopping.\n")
      break

    }#THEN

    state[[loop.counter]] = mb

    if (ia.detect.infinite.loop(mb, state, loop.counter - 1, debug)) {

      # remove the node added to the markov blanket in the last iteration, and
      # do not consider it again for inclusion.
      duplicated.node = mb[length(mb)]
      culprit = c(culprit, duplicated.node)
      mb = mb[mb != duplicated.node]

      if (debug)
        cat("  ! removing", duplicated.node, "from the nodes to test.\n")

      # reset the state list so that no further errors are raised.
      state[[loop.counter]] = NULL
      # reset the loop counter to match.
      loop.counter = loop.counter - 1

      warning("prevented infinite loop in Markov blanket learning (node '", x, "').")

    }#THEN

    # increment the loop counter.
    loop.counter = loop.counter + 1

  }#REPEAT

  return(mb)

}#INTER.IA.MARKOV.BLANKET

