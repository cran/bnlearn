
inter.incremental.association.optimized = function(x, whitelist, blacklist,
  test, alpha, B, strict, debug = FALSE) {

  nodes = names(x)
  mb2 = mb = list()

  # 1. [Compute Markov Blankets]
  for (node in nodes) {

    backtracking = unlist(sapply(mb, function(x){ node %in% x  }))

    mb[[node]] = inter.ia.markov.blanket(node, data = x, nodes = nodes,
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

}#INTER.INCREMENTAL.ASSOCIATION.OPTIMIZED

inter.incremental.association.cluster = function(x, cluster, whitelist,
  blacklist, test, alpha, B, strict, debug = FALSE) {

  nodes = names(x)

  # 1. [Compute Markov Blankets]
  mb = parLapply(cluster, as.list(nodes), inter.ia.markov.blanket, data = x,
         nodes = nodes, alpha = alpha, B = B, whitelist = whitelist,
         blacklist = blacklist, test = test, debug = debug)
  names(mb) = nodes

  # check markov blankets for consistency.
  mb = bn.recovery(mb, nodes = nodes, strict = strict, mb = TRUE, debug = debug)

  # 2. [Compute Graph Structure]
  mb = parLapply(cluster, as.list(nodes), neighbour, mb = mb, data = x,
         alpha = alpha, B = B, whitelist = whitelist, blacklist = blacklist,
         test = test, debug = debug)
  names(mb) = nodes

  # check neighbourhood sets for consistency.
  mb = bn.recovery(mb, nodes = nodes, strict = strict, debug = debug)

  return(mb)

}#INTER.INCREMENTAL.ASSOCIATION.CLUSTER

inter.incremental.association = function(x, whitelist, blacklist, test,
  alpha, B, strict, debug = FALSE) {

  nodes = names(x)

  # 1. [Compute Markov Blankets]
  mb = lapply(as.list(nodes), inter.ia.markov.blanket, data = x,
         nodes = nodes, alpha = alpha, B = B, whitelist = whitelist,
         blacklist = blacklist, test = test, debug = debug)
  names(mb) = nodes

  # check markov blankets for consistency.
  mb = bn.recovery(mb, nodes = nodes, strict = strict, mb = TRUE, debug = debug)

  # 2. [Compute Graph Structure]
  mb = lapply(as.list(nodes), neighbour, mb = mb, data = x, alpha = alpha,
         B = B, whitelist = whitelist, blacklist = blacklist, test = test,
         debug = debug)
  names(mb) = nodes

  # check neighbourhood sets for consistency.
  mb = bn.recovery(mb, nodes = nodes, strict = strict, debug = debug)

  return(mb)

}#INTER.INCREMENTAL.ASSOCIATION

inter.ia.markov.blanket = function(x, data, nodes, alpha, B, whitelist, blacklist,
  start = character(0), backtracking = NULL, test, debug = FALSE) {

  nodes = nodes[nodes != x]
  culprit = known.good = known.bad = c()
  whitelisted = nodes[sapply(nodes,
          function(y) { is.whitelisted(whitelist, c(x, y), either = TRUE) })]
  mb = start
  loop.counter = 1
  state = vector(5 * length(nodes), mode = "list")

  del.node = function(y, x, test) {

    if (debug)
      cat("  * checking node", y, "for exclusion (shrinking phase).\n")

    a = conditional.test(x, y, mb[mb != y], data = data, test = test, B = B,
          alpha = alpha)

    if (a > alpha) {

      if (debug) {

        cat("    > node", y, "removed from the markov blanket. ( p-value:", a, ")\n")
        cat("    > conditioning subset: '", mb[mb != y], "'\n")

      }#THEN

      # update the markov blanket.
      assign("mb", mb[mb != y], envir = sys.frame(-3))

      return(NULL)

    }#THEN
    else if (debug) {

      cat("    > node", y, "remains in the markov blanket. ( p-value:", a, ")\n")

    }#THEN

  }#DEL.NODE

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
  nodes = nodes[!(nodes %in% mb)]
  # blacklist is not checked, not all nodes in a markov blanket must be
  # neighbours.

  # use backtracking for a further screening of the nodes to be checked.
  if (!is.null(backtracking)) {

    # nodes whose markov blanket includes this node are included, because
    # X \in MB(Y) <=> Y \in MB(X)
    known.good = names(backtracking[backtracking])
    mb = unique(c(mb, known.good))

    # and vice versa X \not\in MB(Y) <=> Y \not\in MB(X)
    known.bad = names(backtracking[!backtracking])

    # both are not to be checked for inclusion/exclusion.
    nodes = nodes[!(nodes %in% names(backtracking))]

    if (debug) {

      cat("    * known good (backtracking): '", known.good, "'.\n")
      cat("    * known bad (backtracking): '", known.bad, "'.\n")
      cat("    * nodes still to be tested for inclusion: '", nodes, "'.\n")

    }#THEN

  }#THEN

  repeat {

    # get a snapshot of the markov blanket status.
    mb.snapshot = mb

    # get an association measure for each of the available nodes.
    association = sapply(nodes[!(nodes %in% c(mb, culprit))], conditional.test, x, sx = mb,
                    test = test, data = data, B = B, alpha = alpha)

    # stop if there are no candidates for inclusion; the markov blanket
    # would obviously be unchanged.
    if (all(association > alpha) || length(nodes) == 0 || is.null(nodes)) break

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

    # whitelisted nodes are neighbours, they cannot be removed from the
    # markov blanket; the last node added in phase I will never be removed,
    # because the tests for inclusion and removal are identical.
    # known.good nodes from backtracking are not to be removed, either.
    if (length(mb) > 1)
      sapply(mb[!(mb %in% c(known.good, to.add, whitelisted))], del.node, x = x, test = test)

    if (identical(mb, mb.snapshot)) {

      if (debug)
        cat("  @ markov blanket is unchanged, stopping.")
      break

    }#THEN

    state[[loop.counter]] = mb

    # the stepwise learning technique used here is stateless; make it stateful
    # to detect and break infinite loops.
    duplicated.check = duplicated(state[seq(loop.counter)])

    if (any(duplicated.check)) {

      # remove the node added to the markov blanket in the last iteration, and
      # do not consider it again for inclusion.
      duplicated.node = mb[length(mb)]
      culprit = c(culprit, duplicated.node)
      mb = mb[mb != duplicated.node]

      if (debug) {

        cat("  ! recurring markov blanket configuration detected (",
          state[[loop.counter]], ").\n")
        cat("  ! retracing the steps of the learning process:\n")
        sapply(state[seq(loop.counter)],
          function(str) {
            cat("    >", paste(str, collapse = " "), "\n")
          })
        cat("  ! removing", duplicated.node, "from the nodes to test.\n")

      }#THEN

      # reset the state list so that no further errors are raised.
      state[[loop.counter]] = NULL
      # reset the loop counter to match.
      loop.counter = loop.counter - 1

      warning("prevented infinite loop in Markov blanket learning (node '", x, "').")

    }#THEN

    # increment the loop counter.
    loop.counter = loop.counter + 1

  }#REPEAT

  mb

}#INTER.IA.MARKOV.BLANKET

