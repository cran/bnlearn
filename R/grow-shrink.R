
grow.shrink.optimized = function(x, whitelist, blacklist, test, alpha,
  strict, debug) {

  nodes = names(x)
  mb2 = mb = list()

  # 1. [Compute Markov Blankets]
  for (node in nodes) {

    backtracking = unlist(sapply(mb, function(x){ node %in% x  }))

    mb[[node]] = gs.markov.blanket(node, data = x, nodes = nodes,
         alpha = alpha, whitelist = whitelist, blacklist = blacklist,
         backtracking = backtracking, test = test, debug = debug)

  }#FOR

  # check markov blankets' consistency.
  mb = mb.recovery(mb, nodes = nodes, strict = strict, debug = debug)

  # 2. [Compute Graph Structure]
  for (node in nodes) {

    backtracking = unlist(sapply(mb2, function(x){ node %in% x[["nbr"]]  }))

    # save results in a copy of mb;
    mb2[[node]] = neighbour(node, mb = mb, data = x, alpha = alpha,
         whitelist = whitelist, blacklist = blacklist,
         backtracking = backtracking, test = test, debug = debug)

  }#FOR

  # update mb with the results of neighbour().
  mb = mb2

  # hope it's never called ...
  mb = nbr.recovery(mb, nodes = nodes, strict = strict, debug = debug)

  return(mb)

}#GROW.SHRINK.OPTIMIZED

grow.shrink.cluster = function(x, cluster, whitelist, blacklist, test,
  alpha, strict, debug) {

  nodes = names(x)

  # 1. [Compute Markov Blankets]
  mb = parLapply(cluster, as.list(nodes), gs.markov.blanket, data = x,
         nodes = nodes, alpha = alpha, whitelist = whitelist,
         blacklist = blacklist, test = test, debug = debug)
  names(mb) = nodes

  # check markov blankets' consistency.
  mb = mb.recovery(mb, nodes = nodes, strict = strict, debug = debug)

  # 2. [Compute Graph Structure]
  mb = parLapply(cluster, as.list(nodes), neighbour, mb = mb, data = x,
         alpha = alpha, whitelist = whitelist, blacklist = blacklist,
         test = test, debug = debug)
  names(mb) = nodes

  # hope it's never called ...
  mb = nbr.recovery(mb, nodes = nodes, strict = strict, debug = debug)

  return(mb)

}#GROW.SHRINK.CLUSTER

grow.shrink = function(x, whitelist, blacklist, test, alpha,
  strict, debug) {

  nodes = names(x)

  # 1. [Compute Markov Blankets]
  mb = lapply(as.list(nodes), gs.markov.blanket, data = x, nodes = nodes,
         alpha = alpha, whitelist = whitelist, blacklist = blacklist,
         test = test, debug = debug)
  names(mb) = nodes

  # check markov blankets' consistency.
  mb = mb.recovery(mb, nodes = nodes, strict = strict, debug = debug)

  # 2. [Compute Graph Structure]
  mb = lapply(as.list(nodes), neighbour, mb = mb, data = x, alpha = alpha,
         whitelist = whitelist, blacklist = blacklist, test = test,
         debug = debug)
  names(mb) = nodes

  # hope it's never called ...
  mb = nbr.recovery(mb, nodes = nodes, strict = strict, debug = debug)

  return(mb)

}#GROW.SHRINK

gs.markov.blanket = function(x, data, nodes, alpha, whitelist, blacklist,
  backtracking = NULL, test, debug) {

  nodes = nodes[nodes != x]
  known.good = known.bad = c()
  whitelisted = nodes[sapply(nodes,
          function(y) { is.whitelisted(whitelist, c(x,y), either = TRUE) })]
  mb = c()

  # growing phase
  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* detecting markov blanket of", x, ".\n")

  }#THEN

  # whitelisted nodes are included by default (if there's a direct arc
  # between them of course they are in each other's markov blanket).
  # arc direction is irrelevant here.
  mb = whitelisted
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

  add.node = function(y, x, test) {

    if (debug)
      cat("  * checking node", y, "for inclusion.\n")

    a = conditional.test(x, y, mb, data = data, test = test)

    if (a <= alpha) {

      if (debug) {

        cat("    > node", y, "included in the markov blanket ( p-value:", a, ").\n")
        cat("    > markov blanket now is '", c(mb, y), "'.\n")

      }#THEN

      assign('mb', c(mb, y), envir = sys.frame(-3))

    }#THEN
    else if (debug) {

      cat("    >", x, "indep.", y, "given '", mb, "' ( p-value:", a, ").\n")

    }#THEN

  }#ADD.NODE

  # first pass detects parents and children
  sapply(nodes, add.node, x = x, test = test)
  # second pass detects parents of the son
  remaining = nodes[!(nodes %in% mb)]
  # no second pass is needed if the only node standing is the last one,
  # or if no node was added in the first pass.
  if ((!identical(remaining, nodes)) && (!identical(remaining, nodes[length(nodes)])))
      sapply(remaining, add.node, x = x, test = test)
  # remember the last node added to the markov blanket, which is not to be
  # tested for removal.
  last.added = mb[length(mb)]

  # shrinking phase
  del.node = function(y, x, test) {

    if (debug)
      cat("  * checking node", y, "for exclusion (shrinking phase).\n")

    a = conditional.test(x, y, mb[mb != y], data = data, test = test)

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

  # whitelisted nodes are neighbours, they cannot be removed from the
  # markov blanket; the last node added in phase I will never be removed,
  # because the tests for inclusion and removal are identical.
  # known.good nodes from backtracking are not to be removed, either.
  if (length(mb) > 1)
    sapply(mb[!(mb %in% c(known.good, last.added, whitelisted))], del.node, x = x, test = test)

  mb

}#GS.MARKOV.BLANKET

