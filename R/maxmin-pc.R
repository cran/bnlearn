
maxmin.pc.optimized = function(x, whitelist, blacklist, test, alpha, B,
    max.sx = ncol(x), strict, complete, debug = FALSE, noise.levels = NULL) {

  nodes = names(x)
  mb = list()

  for (node in nodes) {

    backtracking = unlist(sapply(mb, function(x){ node %in% x$nbr }))

    # 1. [Forward Phase (I)]
    mb[[node]] = maxmin.pc.forward.phase(node, data = x, nodes = nodes,
         alpha = alpha, B = B, whitelist = whitelist, blacklist = blacklist,
         backtracking = backtracking, test = test, max.sx = max.sx,
         optimized = TRUE, complete = complete, debug = debug, noise.levels = noise.levels)

    # 2. [Backward Phase (II)]
    mb[[node]] = neighbour(node, mb = mb, data = x, alpha = alpha,
         B = B, whitelist = whitelist, blacklist = blacklist,
         backtracking = backtracking, test = test, max.sx = max.sx,
         markov = FALSE, complete = complete, debug = debug, noise.levels = noise.levels)

  }#FOR

  # make up a set of believable Markov blankets, using all the nodes within
  # distance 2 from the target node (which is a superset).
  for (node in nodes)
    mb[[node]]$mb = fake.markov.blanket(mb, node)

  # check neighbourhood sets for consistency.
  mb = bn.recovery(mb, nodes = nodes, strict = strict, debug = debug)

  return(mb)

}#MAXMIN.PC.OPTIMIZED

maxmin.pc = function(x, cluster = NULL, whitelist, blacklist, test, alpha, B,
  strict, max.sx = ncol(x), complete, debug = FALSE, noise.levels = NULL) {

  nodes = names(x)

  # 1. [Forward Phase (I)]
  mb = smartSapply(cluster, as.list(nodes), maxmin.pc.forward.phase, data = x,
         nodes = nodes, alpha = alpha, B = B, whitelist = whitelist,
         blacklist = blacklist, test = test, max.sx = max.sx,
         optimized = FALSE, complete = complete, debug = debug, noise.levels = noise.levels)
  names(mb) = nodes

  # 2. [Backward Phase (II)]
  mb = smartSapply(cluster, as.list(nodes), neighbour, mb = mb, data = x,
         alpha = alpha, B = B, whitelist = whitelist, blacklist = blacklist,
         test = test, max.sx = max.sx, markov = FALSE, complete = complete,
         debug = debug, noise.levels = noise.levels)
  names(mb) = nodes

  # make up a set of believable Markov blankets, using all the nodes within
  # distance 2 from the target node (which is a superset).
  for (node in nodes)
    mb[[node]]$mb = fake.markov.blanket(mb, node)

  # check neighbourhood sets for consistency.
  mb = bn.recovery(mb, nodes = nodes, strict = strict, debug = debug)

  return(mb)

}#MAXMIN.PC

maxmin.pc.forward.phase = function(x, data, nodes, alpha, B, whitelist,
  blacklist, backtracking = NULL, test, max.sx = ncol(x), optimized = TRUE,
  complete, debug = FALSE, noise.levels = NULL) {

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
  if (!is.null(backtracking) && optimized) {

    # nodes whose neighbourhood includes this node are included, because
    # X \in NBR(Y) <=> Y \in NBR(X); they can be removed in the backward
    # phase later on if they are false positives.
    known.good = names(backtracking[backtracking])
    cpc = unique(c(cpc, known.good))

    # and vice versa X not adiacent to Y <=> Y not adiacent to X
    known.bad = names(backtracking[!backtracking])

    # from empirical observations, MMPC has very few false negatives and
    # therefore it is fairly safe not to test known.bad nodes.
    nodes  = nodes[nodes %!in% known.bad]

    if (debug) {

      cat("    * known good (backtracking): '", known.good, "'.\n")
      cat("    * known bad (backtracking): '", known.bad, "'.\n")
      cat("    * nodes still to be tested for inclusion: '",
        nodes[nodes %!in% cpc], "'.\n")

    }#THEN

  }#THEN

  # phase I (stepwise forward selection)
  repeat {

    # stop testing if the conditioning set grows too large.
    if (length(cpc) > max.sx)
      break

    # do not check nodes which have a p-value above the alpha threshold, as
    # it can only increase; also do not check 'known bad' ones.
    to.be.checked = setdiff(names(which(association < alpha)), c(cpc, known.bad))

    # get an association measure for each of the available nodes.
    association = sapply(to.be.checked, maxmin.pc.heuristic.optimized, y = x,
                    sx = cpc, data = data, test = test, alpha = alpha, B = B,
                    association = association, complete = complete, debug = debug,
                    noise.levels = noise.levels)

    # stop if there are no candidates for inclusion.
    if (all(association > alpha) || length(nodes) == 0 || is.null(nodes)) break
    # get the one which maximizes the association measure.
    to.add = names(which.min(association))

    if (debug) {

      cat("  @", to.add, "accepted as a parent/children candidate ( p-value:",
        association[to.add], ").\n")
      cat("  > current candidates are '", c(cpc, to.add), "'.\n")

    }#THEN

    if (association[to.add] <= alpha) {

      cpc = c(cpc, to.add)
      nodes = nodes[nodes != to.add]

    }#THEN

  }#REPEAT

  return(cpc)

}#MAXMIN.PC.FORWARD.PHASE

maxmin.pc.heuristic.optimized = function(x, y, sx, data, test, alpha, B,
    association, complete, debug = FALSE, noise.levels = NULL) {

  min.assoc = association[x]

  if (debug) {

    cat("  * checking node", x ,"for association.\n")
    cat("    > starting with association", min.assoc, ".\n")

  }#THEN

  # generate only the subsets of the current parent/children set which include
  # node added last; the rest are considered to be already tested against.
  last = sx[length(sx)]
  sx = sx[-length(sx)]
  
  new.min.assoc = allsubs.test(x = x, y = y, sx = sx, fixed = last, data = data,
                    test = test, B = B, alpha = alpha, complete = complete,
                    debug = debug, noise.levels = noise.levels)

  min.assoc = max(min.assoc, new.min.assoc["max.p.value"])

  if (debug)
    cat("    > node", x, "has a minimum association of", min.assoc, ".\n")

  return(min.assoc)

}#MAXMIN.PC.HEURISTIC.OPTIMIZED

