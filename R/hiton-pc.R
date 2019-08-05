
si.hiton.pc.backend = function(x, cluster = NULL, whitelist, blacklist,
  test, alpha, B, max.sx = ncol(x), complete, debug = FALSE) {

  nodes = names(x)

  # 1. [Forward Phase (I)]
  mb = smartSapply(cluster, as.list(nodes), si.hiton.pc.heuristic, data = x,
         nodes = nodes, alpha = alpha, B = B, whitelist = whitelist,
         blacklist = blacklist, test = test, max.sx = max.sx,
         complete = complete, debug = debug)
  names(mb) = nodes

  # 2. [Backward Phase (II)]
  mb = smartSapply(cluster, as.list(nodes), neighbour, mb = mb, data = x,
         alpha = alpha, B = B, whitelist = whitelist, blacklist = blacklist,
         test = test, debug = debug, max.sx = max.sx, empty.dsep = FALSE,
         markov = FALSE, complete = complete)
  names(mb) = nodes

  # make up a set of believable Markov blankets, using all the nodes within
  # distance 2 from the target node (which is a superset).
  for (node in nodes)
    mb[[node]]$mb = fake.markov.blanket(mb, node)

  # check neighbourhood sets for consistency.
  mb = bn.recovery(mb, nodes = nodes, debug = debug)

  return(mb)

}#SI.HITON.PC.BACKEND

si.hiton.pc.heuristic = function(x, data, nodes, alpha, B, whitelist, blacklist,
    test, max.sx = ncol(x), complete, debug = FALSE) {

  nodes = nodes[nodes != x]
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

   # no nodes to check, nothing to do, move along.
  if (length(nodes) == 0)
    return(cpc)

  # get a marginal association measure for each of the available nodes.
  association = indep.test(nodes, x, sx = character(0), test = test,
                  data = data, B = B, alpha = alpha, complete = complete)

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

  # phase I (stepwise forward selection).
  repeat {

    # stop if there are no candidates for inclusion.
    if (all(association > alpha) || length(nodes) == 0 || is.null(nodes) ||
          length(cpc) > max.sx)
      break
    # get the one which maximizes the association measure.
    to.add = names(which.min(association))

    # check whether the node is independent of the target given a subset of
    # the current Markov blanket.
    candidate = si.hiton.pc.backward(target = x, candidate = to.add, cpc = cpc,
                  data = data, test = test, alpha = alpha, B = B,
                  complete = complete, debug = debug)

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
si.hiton.pc.backward = function(target, candidate, cpc, data, test, alpha, B,
    complete, debug) {

  # the nodes are always marginally associated, otherwise the candidate would
  # not have been chosen as such.
  if (length(cpc) == 0)
    return(TRUE)

  if (debug)
    cat("* backward phase for candidate node", candidate, ".\n")

  a = allsubs.test(x = target, y = candidate, sx = cpc, min = 1L, data = data,
        test = test, alpha = alpha, B = B, complete = complete, debug = debug)

  return(a["p.value"] <= alpha)

}#SI.HITON.PC.BACKWARD

