
fast.incremental.association.optimized = function(x, whitelist, blacklist,
  test, alpha, B, strict, debug = FALSE) {

  nodes = names(x)
  mb2 = mb = list()

  # 1. [Compute Markov Blankets]
  for (node in nodes) {

    backtracking = unlist(sapply(mb, function(x){ node %in% x  }))

    mb[[node]] = fast.ia.markov.blanket(node, data = x, nodes = nodes,
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

}#FAST.INCREMENTAL.ASSOCIATION.OPTIMIZED

fast.incremental.association = function(x, cluster = NULL, whitelist,
  blacklist, test, alpha, B, strict, debug = FALSE) {

  nodes = names(x)

  # 1. [Compute Markov Blankets]
  mb = smartLapply(cluster, as.list(nodes), fast.ia.markov.blanket,
         data = x, nodes = nodes, alpha = alpha, B = B, whitelist = whitelist,
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

}#FAST.INCREMENTAL.ASSOCIATION

fast.ia.markov.blanket = function(x, data, nodes, alpha, B, whitelist, blacklist,
  start = character(0), backtracking = NULL, test, debug = FALSE) {

  nodes = nodes[nodes != x]
  known.good = known.bad = c()
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

  # use backtracking for a further screening of the nodes to be checked.
  if (!is.null(backtracking)) {

    # nodes whose markov blanket includes this node are included, because
    # X \in MB(Y) <=> Y \in MB(X); they can be removed in the backward phase
    # later on if they are false positives.
    known.good = names(backtracking[backtracking])
    mb = unique(c(mb, known.good))

    # and vice versa X \not\in MB(Y) <=> Y \not\in MB(X)
    known.bad = names(backtracking[!backtracking])

    if (debug) {

      cat("    * known good (backtracking): '", known.good, "'.\n")
      cat("    * known bad (backtracking): '", known.bad, "'.\n")
      cat("    * nodes still to be tested for inclusion: '",
        nodes[nodes %!in% mb], "'.\n")

    }#THEN

  }#THEN

  repeat {

    # stop if there are no nodes left.
    if (length(nodes[nodes %!in% mb]) == 0 || is.null(nodes))
      break

    # growing phase.
    # reset the insufficient.data boolean flag.
    insufficient.data = FALSE

    # get an association measure for each of the available nodes.
    association = indep.test(nodes[nodes %!in% mb], x, sx = mb, test = test,
                    data = data, B = B, alpha = alpha)

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
    # markov blanket; on the other hand, known.good nodes from backtracking
    fixed = whitelisted[whitelisted != ""]

    pv = roundrobin.test(x = x, z = mb, fixed = fixed, data = data, test = test,
           B = B, alpha = alpha, debug = debug)

    mb = intersect(mb, c(names(pv[pv < alpha]), fixed))

    # if there are not enough observations and no new node has been included
    # in the markov blanket, stop iterating.
    if (insufficient.data && (mb.old.length == length(mb))) break

    # do not touch the check the nodes in the markov blanket again.
    nodes = nodes[nodes %!in% mb]

  }#REPEAT

  return(mb)

}#FAST.IA.MARKOV.BLANKET

