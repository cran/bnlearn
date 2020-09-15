
hybrid.pc.backend = function(x, cluster = NULL, whitelist, blacklist,
  test, alpha, B, max.sx = ncol(x), complete, debug = FALSE) {

  nodes = names(x)

  mb = smartSapply(cluster, as.list(nodes), hybrid.pc.heuristic, data = x,
         nodes = nodes, alpha = alpha, B = B, whitelist = whitelist,
         blacklist = blacklist, test = test, max.sx = max.sx,
         complete = complete, debug = debug)
  names(mb) = nodes

  # make up a set of believable Markov blankets, using all the nodes within
  # distance 2 from the target node (which is a superset).
  for (node in nodes)
    mb[[node]]$mb = fake.markov.blanket(mb, node)

  # check neighbourhood sets for consistency.
  mb = bn.recovery(mb, nodes = nodes, debug = debug)

  return(mb)

}#HYBRID.PC.BACKEND

hybrid.pc.heuristic = function(x, data, nodes, alpha, B, whitelist, blacklist,
    test, max.sx = ncol(data), complete, debug = FALSE) {

  # identify the parents-and-children superset.
  pvalues = hybrid.pc.de.pcs(x = x, data = data, nodes = nodes, alpha = alpha,
              B = B, whitelist = whitelist, blacklist = blacklist,
              test = test, complete = complete, debug = debug)
  pc.superset = names(pvalues)

  # if the superset contains zero or just one nodes, there is nothing else to do
  # (the superset is not super, it is the right set).
  if (length(pc.superset) < 2)
    return(list(nbr = pc.superset, mb = NULL))

  # identify the spouses superset (some spouses may already be in the
  # parents-and-children superset).
  sp.superset = hybrid.pc.de.sps(x = x, data = data, nodes = nodes,
                  pc.superset = pc.superset, dsep.set = attr(pvalues, "dsep.set"),
                  alpha = alpha, B = B, test = test, max.sx = max.sx,
                  complete = complete, debug = debug)

  # if there are just two nodes in the parents-and-children set and no spouse,
  # the superset is necessarily the same as the set..
  if ((length(pc.superset) == 2) && (length(sp.superset) == 0))
    return(list(nbr = pc.superset, mb = pc.superset))

  # the two nodes with the smallest p-values would be examined again using the
  # same low-order tests as before, so just include them.
  start = names(sort(pvalues))[1:min(length(pvalues), 2)]

  # identify the real parents and children from the supersets.
  pc = hybrid.pc.nbr.search(x, data = data,
         nodes = c(x, pc.superset, sp.superset), alpha = alpha, B = B,
         whitelist = whitelist, blacklist = blacklist, test = test,
         max.sx = max.sx, complete = complete, start = start, debug = debug)

  # one more scan to identify possible false negatives.
  for (node in setdiff(pc.superset, pc)) {

    fn = hybrid.pc.nbr.search(node, data = data,
            nodes = c(x, pc.superset, sp.superset), alpha = alpha, B = B,
            whitelist = whitelist, blacklist = blacklist, test = test,
            max.sx = max.sx, complete = complete, start = start,
            debug = debug, looking.for = x)

    # add the nodes which appear to be neighbours.
    if (x %in% fn) {

      pc = c(pc, node)
      mb = c(mb, node)

      if (debug)
        cat("  @", node, "added to the parents and children. (HPC's OR)\n")

    }#THEN

  }#FOR

  return(list(nbr = pc, mb = c(pc.superset, sp.superset)))

}#HYBRID.PC.HEURISTIC

hybrid.pc.nbr.search = function(x, data, nodes, alpha, B, whitelist, blacklist,
    test, max.sx = ncol(data), complete, debug = FALSE, start = start,
    looking.for = NULL) {

  mb = ia.fdr.markov.blanket(x, data = data, nodes = nodes, alpha = alpha,
         B = B, whitelist = whitelist, blacklist = blacklist, start = start,
         test = test, max.sx = max.sx, complete = complete, debug = debug)

  # if the node is not in the markov blanket it cannot be a neighbour either.
  if (!is.null(looking.for) && (looking.for %!in% mb))
    return(NULL)

  pc = hybrid.pc.filter(x, pc.superset = mb, sp.superset = NULL, data = data,
         alpha = alpha, B = B, whitelist = whitelist, blacklist = blacklist,
         test = test, max.sx = max.sx, complete = complete, debug = debug)

  return(pc)

}#HYBRID.PC.NBR.SEARCH

hybrid.pc.de.pcs = function(x, data, nodes, alpha, B, whitelist, blacklist,
    test, complete, debug = FALSE) {

  whitelisted = nodes[sapply(nodes,
          function(y) { is.whitelisted(whitelist, c(x, y), either = TRUE) })]
  blacklisted = nodes[sapply(nodes,
          function(y) { is.blacklisted(blacklist, c(x, y), both = TRUE) })]

  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* learning the parents and children superset of", x, ".\n")

  }#THEN

  if (debug)
    cat("  * nodes to be tested for inclusion: '",
        nodes[nodes %!in% x], "'.\n")

  # all nodes are candidates initially, except for those whose status is already
  # determined (from whitelist and blacklist).
  to.check = setdiff(nodes, c(x, whitelisted, blacklisted))

  # exclude nodes that are marginally independent from the target.
  association = indep.test(to.check, x, sx = character(0), data = data,
                  test = test, B = B, alpha = alpha, complete = complete)

  to.keep = names(association[association <= alpha])
  to.drop = names(association[association > alpha])
  pvalues = association[to.keep]

  if (debug) {

    cat("  * checking nodes for association.\n")

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

  # sort the candidates in order of increasing association, so that nodes with
  # weak associations are checked for exclusion first.
  pvalues = sort(pvalues, decreasing = TRUE)
  to.check = names(pvalues)
  fixed = whitelisted

  if (debug)
    cat("  * nodes to be tested for exclusion: '", to.check, "'.\n")

  pvalues = structure(c(rep(0, length(fixed)), pvalues),
              names = c(fixed, names(pvalues)))
  dsep.set = list()

  # if there is only a single node left to check, it would be conditional on
  # the empty set which would just repeat an earlier test; nothing left to do.
  if (length(to.check) == 1)
    return(structure(pvalues, dsep.set = dsep.set))

  # exlcude nodes that are independent given a single conditioning node.
  for (node in to.check) {

    # sort the candidates in order of decreasing association, so that nodes with
    # strong associations are checked first.
    to.check.against = setdiff(names(sort(pvalues, decreasing = FALSE)), node)

    if (length(to.check.against) == 0)
      next

    a = allsubs.test(x = x, y = node, sx = to.check.against, min = 1,
          max = 1, data = data, test = test, alpha = alpha, B = B,
          complete = complete, debug = debug)

    if (a["p.value"] > alpha) {

      pvalues = pvalues[names(pvalues) != node]
      dsep.set[[node]] = attr(a, "dsep.set")

    }#THEN
    else {

      pvalues[node] = max(pvalues[node], a["max.p.value"])

    }#ELSE

  }#FOR

  return(structure(pvalues, dsep.set = dsep.set))

}#HYBRID.PC.DE.PCS

hybrid.pc.de.sps = function(x, data, nodes, pc.superset, dsep.set, alpha, B,
    test, max.sx, complete, debug = FALSE) {

  spouses.superset = character(0)

  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* learning the spouses superset of", x, ".\n")
    cat("  * nodes still to be tested for inclusion:",
      nodes[nodes %!in% c(x, pc.superset)], "\n")

  }#THEN

  for (cpc in pc.superset) {

    pvalues = numeric(0)

    # forward selection.
    for (y in setdiff(nodes, c(x, pc.superset))) {

      # if the candidate node d-separates the current node from a node that is
      # not in the superset, it is potentially in the markov blanket and thus a
      # potential spouse.
      if (cpc %in% dsep.set[[y]])
        next

      # skip tests whose conditioning sets are too large (assuming independence
      # means not adding the node to the superset).
      if (length(c(dsep.set[[y]], cpc)) > max.sx)
        next

      if (debug)
        cat("  > checking node", y, "for inclusion.\n")

      a = indep.test(x = x, y = y, sx = c(dsep.set[[y]], cpc), data = data,
            test = test, B = B, alpha = alpha, complete = complete)

      if (debug)
        cat("    > node", x, "is",
          ifelse(a > alpha, "independent from", "dependent on"), y, "given",
          c(dsep.set[[y]], cpc), " ( p-value:", a, ").\n")

      if (a <= alpha) {

        pvalues[y] = a

        if (debug)
          cat("    @ node", y, "added to the spouses superset.\n")

      }#THEN

    }#FOR

    # sort the candidates in order of increasing association, so that nodes with
    # weak associations are checked for exclusion first.
    pvalues = sort(pvalues, decreasing = TRUE)

    # backward selection, to remove false positives.
    for (y in names(pvalues)) {

      sx = setdiff(names(pvalues), y)

      if (length(sx) == 0)
        next

      if (debug)
        cat("  > checking node", y, "for removal.\n")

      a = allsubs.test(x = x, y = y, sx = sx, fixed = cpc, data = data,
            test = test, B = B, alpha = alpha, complete = complete, min = 1,
            max = 1, debug = debug)

      if (a["p.value"] > alpha) {

        pvalues = pvalues[names(pvalues) != y]

        if (debug)
          cat("    @ node", y, "removed from the spouses superset.\n")

      }#THEN

    }#FOR

    spouses.superset = union(spouses.superset, names(pvalues))

  }#FOR

  return(spouses.superset)

}#HYBRID.PC.DE.SPS

hybrid.pc.filter = function(x, pc.superset, sp.superset, data, alpha, B = B,
     whitelist, blacklist, test, max.sx, complete, debug = FALSE) {

  nodes = names(data)
  mb.superset = c(pc.superset, sp.superset)

  whitelisted = nodes[sapply(nodes,
          function(y) { is.whitelisted(whitelist, c(x, y), either = TRUE) })]
  blacklisted = nodes[sapply(nodes,
          function(y) { is.blacklisted(blacklist, c(x, y), both = TRUE) })]

  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* filtering parents and children of", x, ".\n")
    cat("  * blacklisted nodes: '", blacklisted, "'\n")
    cat("  * whitelisted nodes: '", whitelisted, "'\n")
    cat("  * starting with neighbourhood superset: '", pc.superset, "'\n")
    cat("  * with spouses superset: '", sp.superset, "'\n")

  }#THEN

  # make sure blacklisted nodes are not included, and add whitelisted nodes.
  pc.superset = union(setdiff(pc.superset, blacklisted), whitelisted)
  mb.superset = union(mb.superset, whitelisted)

  # if the markov blanket is empty, the neighbourhood is empty as well.
  if (length(mb.superset) == 0)
    return(character(0))

  nbr = function(node) {

    a = allsubs.test(x = x, y = node, sx = setdiff(mb.superset, node),
          data = data, test = test, B = B, alpha = alpha, max = max.sx,
          complete = complete, debug = debug)

    return(as.logical(a["p.value"] < alpha))

  }#NBR

  # do not even try to remove whitelisted nodes.
  pc = names(which(sapply(setdiff(pc.superset, whitelisted), nbr)))

  # make sure whitelisted nodes are always included.
  pc = unique(c(pc, whitelisted))

  return(pc)

}#HYBRID.PC.FILTER

