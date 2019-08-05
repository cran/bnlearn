
fast.incremental.association = function(x, cluster = NULL, whitelist,
  blacklist, test, alpha, B, max.sx = ncol(x), complete, debug = FALSE) {

  nodes = names(x)

  # 1. [Compute Markov Blankets]
  mb = smartSapply(cluster, as.list(nodes), fast.ia.markov.blanket,
         data = x, nodes = nodes, alpha = alpha, B = B, whitelist = whitelist,
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

}#FAST.INCREMENTAL.ASSOCIATION

fast.ia.markov.blanket = function(x, data, nodes, alpha, B, whitelist,
  blacklist, start = character(0), test, max.sx = ncol(x), complete,
  debug = FALSE) {

  nodes = nodes[nodes != x]
  whitelisted = nodes[sapply(nodes,
          function(y) { is.whitelisted(whitelist, c(x, y), either = TRUE) })]
  mb = start
  insufficient.data = FALSE

  # growing phase.
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

    # stop if there are no nodes left, or if we cannot add any more nodes
    # because the conditioning set has grown too large.
    if (length(nodes[nodes %!in% mb]) == 0 || is.null(nodes))
      break

    if (length(mb) > max.sx) {

       if (debug)
         cat("  @ limiting conditioning sets to", max.sx, "nodes.\n")

      break

    }#THEN

    # get a snapshot of the markov blanket status.
    mb.snapshot = mb

    # growing phase.
    # reset the insufficient.data boolean flag.
    insufficient.data = FALSE

    # get an association measure for each of the available nodes.
    association = indep.test(nodes[nodes %!in% mb], x, sx = mb, test = test,
                    data = data, B = B, alpha = alpha, complete = complete)

    if (debug) {

      cat("  * checking nodes for association.\n")
      sapply(names(association),
        function(x) {  cat("    >", x, "has p-value", association[x], ".\n")})

    }#THEN

    # stop if there are no candidates for inclusion.
    if (all(association > alpha))
      break

    # sort the candidates in increasing p-value order.
    association = sort(association[which(association <= alpha)])

    for (node in names(association)) {

      # when speculatively including nodes, avoid forming a markov blanket large
      # enough that the tests in the exclusion phase have conditioning sets
      # larger than the threshold.
      if (length(mb) > max.sx) {

        if (debug)
          cat("    @ skipping node", node,
              "to avoid conditioning tests with conditioning sets larger than",
              max.sx, ".\n")

        next

      }#THEN

      opc = obs.per.cell(x, node, mb, data)

      if ((test %!in% asymptotic.tests) || (opc >= 5)) {

        if (debug) {

          if (test %in% available.continuous.tests) {

            cat("    @", node, "included in the markov blanket (p-value:",
              association[node], ").\n")

          }#THEN
          else {

            cat("    @", node, "included in the markov blanket (p-value:",
              association[node], ", obs/cell:", opc, ").\n")

          }#ELSE
          cat("    > markov blanket (", length(mb) + 1, " nodes ) now is '", c(mb, node), "'.\n")

        }#THEN

        # speculatively add the node if the asymptotic behaviour of the
        # statistical tests is still good.
        mb = c(mb, node)

      }#THEN
      else {

        if (debug)
          cat("  @ not enough observations per cell (", opc ,"), skipping.\n")

        # do not add new nodes if that compromises the asymptotic behaviour
        # of the statistical tests.
        insufficient.data = TRUE
        break

      }#ELSE

    }#FOR

    # shrinking phase.
    mb.old.length = length(mb)

    if (debug)
      cat("  * checking nodes for exclusion.\n")

    # whitelisted nodes are neighbours, they cannot be removed from the
    # markov blanket.
    fixed = whitelisted[whitelisted != ""]

    pv = roundrobin.test(x = x, z = mb, fixed = fixed, data = data, test = test,
           B = B, alpha = alpha, complete = complete, debug = debug)

    mb = intersect(mb, c(names(pv[pv < alpha]), fixed))

    # if the markov blanket did not change, stop iterating.
    if (identical(mb, mb.snapshot)) {

      if (debug)
        cat("  @ markov blanket is unchanged, stopping.\n")
      break

    }#THEN

    # if there are not enough observations and no new node has been included
    # in the markov blanket, stop iterating.
    if (insufficient.data && (mb.old.length == length(mb))) break

    # do not touch the nodes in the markov blanket again.
    nodes = nodes[nodes %!in% mb]

  }#REPEAT

  return(mb)

}#FAST.IA.MARKOV.BLANKET

