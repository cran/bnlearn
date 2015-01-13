
si.hiton.pc.optimized = function(x, whitelist, blacklist, test,
  alpha, B, strict, debug = FALSE) {

  nodes = names(x)
  mb = list()

  for (node in nodes) {

    backtracking = unlist(sapply(mb, function(x){ node %in% x$nbr }))

    # 1. [Forward Phase (I)]
    mb[[node]] = si.hiton.pc.heuristic(node, data = x, nodes = nodes,
         alpha = alpha, B = B, whitelist = whitelist, blacklist = blacklist,
         backtracking = backtracking, test = test, debug = debug)

    # 2. [Backward Phase (II)]
    mb[[node]] = neighbour(node, mb = mb, data = x, alpha = alpha,
         B = B, whitelist = whitelist, blacklist = blacklist,
         backtracking = backtracking, test = test, debug = debug,
         empty.dsep = FALSE, markov = FALSE)

  }#FOR

  # check neighbourhood sets for consistency.
  mb = bn.recovery(mb, nodes = nodes, strict = strict, debug = debug)

  return(mb)

}#SI.HITON.PC.OPTIMIZED

si.hiton.pc.backend = function(x, cluster = NULL, whitelist, blacklist,
  test, alpha, B, strict, debug = FALSE) {

  nodes = names(x)

  # 1. [Forward Phase (I)]
  mb = smartLapply(cluster, as.list(nodes), si.hiton.pc.heuristic, data = x,
         nodes = nodes, alpha = alpha, B = B, whitelist = whitelist,
         blacklist = blacklist, test = test, debug = debug)
  names(mb) = nodes

  # 2. [Backward Phase (II)]
  mb = smartLapply(cluster, as.list(nodes), neighbour, mb = mb, data = x,
         alpha = alpha, B = B, whitelist = whitelist, blacklist = blacklist,
         test = test, debug = debug, empty.dsep = FALSE, markov = FALSE)
  names(mb) = nodes

  # check neighbourhood sets for consistency.
  mb = bn.recovery(mb, nodes = nodes, strict = strict, debug = debug)

  return(mb)

}#SI.HITON.PC.BACKEND

si.hiton.pc.heuristic = function(x, data, nodes, alpha, B, whitelist, blacklist,
    backtracking = NULL, test, debug = FALSE) {

  nodes = nodes[nodes != x]
  known.good = known.bad = c()
  whitelisted = nodes[sapply(nodes,
          function(y) { is.whitelisted(whitelist, c(x, y), either = TRUE) })]
  blacklisted = nodes[sapply(nodes,
          function(y) { is.blacklisted(blacklist, c(x, y), both = TRUE) })]
  cpc = c()
  association = structure(numeric(length(nodes)), names = nodes)
  to.add = ""

  # growing phase
  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* forward phase for node", x, ".\n")

  }#THEN

  # whitelisted nodes are included, and blacklisted nodes are excluded.
  cpc = whitelisted
  nodes = nodes[nodes %!in% c(cpc, blacklisted)]

  # use backtracking for a further screening of the nodes to be checked.
  if (!is.null(backtracking)) {

    # X adiacent to Y <=> Y adiacent to X
    known.good = names(backtracking[backtracking])
    cpc = unique(c(cpc, known.good))

    # and vice versa X not adiacent to Y <=> Y not adiacent to X
    known.bad = names(backtracking[!backtracking])

    # known.good nodes are not to be checked for inclusion, and the "nodes"
    # is resetted below so we can just remove them.
    nodes = nodes[nodes %!in% known.good]

    if (debug) {

      cat("    * known good (backtracking): '", known.good, "'.\n")
      cat("    * known bad (backtracking): '", known.bad, "'.\n")
      cat("    * nodes still to be tested for inclusion: '", nodes, "'.\n")

    }#THEN

    # check whether known.good nodes are false positives by running an ad-hoc
    # backward step.
    for (cpn in known.good) {

      candidate = si.hiton.pc.backward(target = x, candidate = cpn,
                    cpc = cpc[cpc != cpn], data = data, test = test,
                    alpha = alpha, B = B, debug = debug)

      if (candidate) {

        if (debug) {

          cat("  @", cpn, "accepted as a parent/children candidate.\n")
          cat("  > current candidates are '", cpc, "'.\n")

        }#THEN

      }#THEN
      else {

        # drop this node, it's apparently a false positive.
        cpc = cpc[cpc != cpn]

        if (debug) {

          cat("  @", cpn, "rejected as a parent/children candidate.\n")
          cat("  > current candidates are '", cpc, "'.\n")

        }#THEN

      }#ELSE

    }#FOR

  }#THEN

  # no nodes to check, nothing to do, move along.
  if (length(nodes) == 0)
    return(cpc)

  # get a marginal association measure for each of the available nodes.
  association = indep.test(nodes, x, sx = character(0), test = test,
                  data = data, B = B, alpha = alpha)

  to.keep = names(association[association <= alpha])
  to.drop = names(association[association > alpha])

  if (debug) {

    cat("  * checking nodes for association.\n")
    cat("  > starting with neighbourhood '", cpc, "'.\n")

    if (length(to.keep) > 0) {

      cat("  * nodes that are still candidates for inclusion.\n")
      sapply(to.keep,
        function(x) {  cat("    >", x, "has p-value", association[x], ".\n")})

    }#THEN

    if (length(to.drop) > 0) {

      cat("  * nodes that will be disregarded from now on.\n")
      sapply(to.drop,
        function(x) {  cat("    >", x, "has p-value", association[x], ".\n")})

    }#THEN

  }#THEN

  # keep around only the nodes that have a significant marginal association.
  nodes = nodes[nodes %in% names(association[association <= alpha])]

  # stop if there are no candidates for inclusion.
  if (all(association > alpha))
    return(cpc)

  # phase I (stepwise forward selection)
  repeat {

    # stop if there are no candidates for inclusion.
    if (all(association > alpha) || length(nodes) == 0 || is.null(nodes)) break
    # get the one which maximizes the association measure.
    to.add = names(which.min(association))

    # check whether the node is independent of the target given a subset of
    # the current Markov blanket.
    candidate = si.hiton.pc.backward(target = x, candidate = to.add, cpc = cpc,
                  data = data, test = test, alpha = alpha, B = B, debug = debug)

    if (candidate) {

      if (debug) {

        cat("  @", to.add, "accepted as a parent/children candidate ( p-value:",
          association[to.add], ").\n")
        cat("  > current candidates are '", c(cpc, to.add), "'.\n")

      }#THEN

      # add the node to the candidate parents-children set.
      cpc = c(cpc, to.add)

    }#THEN

    # remove it from the set of the nodes under consideration.
    nodes = nodes[nodes != to.add]
    association = association[names(association) != to.add]

  }#REPEAT

  return(cpc)

}#SI.HITON.PC.HEURISTIC

# backward stage of HITON-PC.
si.hiton.pc.backward = function(target, candidate, cpc, data, test, alpha, B, debug) {

  # the nodes are always marginally associated, otherwise the candidate would not
  # have been chosen as such.
  if (length(cpc) == 0)
    return(TRUE)

  if (debug)
    cat("* backward phase for candidate node", candidate, ".\n")

  allsubs.test(x = target, y = candidate, sx = cpc, min = 1L, data = data,
    test = test, alpha = alpha, B = B, debug = debug)[1] <= alpha

}#SI.HITON.PC.BACKWARD

