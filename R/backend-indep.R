
# common steps in constraint-based structure learning: infer arc orientation
# from the skeleton and the data.
learn.arc.directions = function(x, cluster = NULL, local.structure, whitelist,
    blacklist, test, alpha, B = NULL, data, max.sx = ncol(data), complete,
    debug = FALSE) {

  nodes = names(x)

  # build a list of the undirected arcs in the graph.
  arcs = mb2arcs(local.structure, nodes)

  # apply blacklist to the arc set.
  to.drop = !apply(arcs, 1, function(x){ is.blacklisted(blacklist, x) })
  arcs = arcs[to.drop, , drop = FALSE]

  # 3. [Orient Edges]
  # 3.1 detect v-structures.
  vs = do.call("rbind",
         vstruct.detect(nodes = nodes, arcs = arcs, mb = local.structure,
           data = x, alpha = alpha, B = B, test = test, blacklist = blacklist,
           max.sx = max.sx, complete = complete, debug = debug))
  rownames(vs) = NULL

  if (!is.null(vs)) {

    # 3.2 sort them in p-value order.
    vs = vs[order(vs[,"max_a"], decreasing = FALSE),]
    # 3.3 apply.
    arcs = vstruct.apply(arcs = arcs, vs = vs, nodes = nodes, debug = debug)

  }#THEN

  # save the status of the learning algorithm.
  learning = list(whitelist = whitelist, blacklist = blacklist,
    test = test, args = list(alpha = alpha), ntests = test.counter())
  # include also the number of permutations/bootstrap samples if it makes sense.
  if (!is.null(B))
    learning$args$B = B

  # 4. propagate directions.
  pdag = list(learning = learning,
              nodes = structure(rep(0, length(nodes)), names = nodes),
              arcs = arcs)

  return(cpdag.backend(pdag, moral = TRUE, fix = TRUE, debug = debug))

}#SECOND.PRINCIPLE

# construct a fake markov blanket using all the nodes within distance 2.
fake.markov.blanket = function(learn, target) {

  mb = c(unlist(lapply(learn[[target]]$nbr,
         function(cur) learn[[cur]]$nbr)), learn[[target]]$nbr)
  mb = setdiff(unique(mb), target)

  return(mb)

}#FAKE.MARKOV.BLANKET

# build the neighbourhood of a node from the markov blanket.
neighbour = function(x, mb, data, alpha, B = NULL, whitelist, blacklist,
  test, empty.dsep = TRUE, markov = TRUE, max.sx = ncol(x), complete,
  debug = FALSE) {

  # initialize the neighbourhood using the markov blanket.
  candidate.neighbours = mb[[x]]

  # if the markov blanket is empty there's nothing to do.
  if (length(candidate.neighbours) == 0) {

    if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("* markov blanket of", x, "is empty; the neighbourhood as well.\n")

    }#THEN

    return(list(mb = character(0), nbr = character(0)))

  }#THEN

  # whitelisted nodes are included (arc orientation is irrelevant), and
  # blacklisted nodes are removed if both directed arcs are banned and both
  # are not in the whitelist.
  blacklisted = candidate.neighbours[sapply(candidate.neighbours,
          function(y) { is.blacklisted(blacklist, c(x, y), both = TRUE) })]
  whitelisted = candidate.neighbours[sapply(candidate.neighbours,
          function(y) { is.whitelisted(whitelist, c(x, y), either = TRUE) })]

  candidate.neighbours = setdiff(candidate.neighbours, blacklisted)
  candidate.neighbours = union(candidate.neighbours, whitelisted)

  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* learning neighbourhood of", x, ".\n")
    cat("  * blacklisted nodes: '", blacklisted, "'\n")
    cat("  * whitelisted nodes: '", whitelisted, "'\n")
    cat("  * starting with neighbourhood: '", candidate.neighbours, "'\n")

  }#THEN

  # nothing much to do, just return.
  if (length(candidate.neighbours) <= 1)
    return(list(mb = mb[[x]], nbr = candidate.neighbours))

  # do not even try to remove whitelisted nodes; on the other hand.
  for (y in setdiff(candidate.neighbours, whitelisted)) {

    if (debug)
      cat("  * checking node", y, "for neighbourhood.\n")

    # choose the smaller set of possible d-separating nodes.
    if (markov)
      dsep.set = smaller(setdiff(mb[[x]], y), setdiff(mb[[y]], x))
    else
      dsep.set = setdiff(mb[[x]], y)

    if (debug)
      cat("    > dsep.set = '", dsep.set, "'\n")

    a = allsubs.test(x = x, y = y, sx = dsep.set,
          min = ifelse(empty.dsep, 0, 1), max = min(length(dsep.set), max.sx),
          data = data, test = test, alpha = alpha, B = B,
          complete = complete, debug = debug)

    # update the neighbourhood.
    if (a["p.value"] > alpha)
      candidate.neighbours = candidate.neighbours[candidate.neighbours != y]

    if (debug)
      cat("    > node", y, "is", ifelse(a["p.value"] > alpha, "not", "still"),
        "a neighbour of", x, ". ( p-value:", a["p.value"], ")\n")

  }#FOR

  return(list(mb = mb[[x]], nbr = candidate.neighbours))

}#NEIGHBOUR

# detect v-structures in the graph.
vstruct.detect = function(nodes, arcs, mb, data, alpha, B = NULL, test,
    blacklist, max.sx = ncol(data), complete, debug = FALSE) {

  vstruct.centered.on = function(x, mb, data, dsep.set) {

    if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("* v-structures centered on", x, ".\n")

    }#THEN

    tos = arcs[(arcs[, "to"] == x), "from"]

    if (length(tos) < 2)
      return(NULL)

    # build a list of possibile parents for the node x, i.e. all the subsets
    # of size 2 of the nodes connected to x by incoming arcs.
    tos.combs = subsets(tos, 2)
    vs = NULL

    for (j in 1:nrow(tos.combs)) {

      y = tos.combs[j, 1]
      z = tos.combs[j, 2]

      if (debug)
        cat("  * checking", y, "->", x, "<-", z, "\n")

      # check there's no arc from y to z and vice versa.
      if (is.listed(arcs, c(y, z), either = TRUE))
        next
      # check whether this is because the arc is blacklisted; d-separating sets
      # are not learned at all in that case.
      if (is.null(blacklist))
        blacklisted.arc = FALSE
      else
        blacklisted.arc = is.listed(blacklist, c(y, z), both = TRUE)

      # if d-separating sets have been saved during the first part of structure
      # learning, it is possible to detect v-structures without additional tests;
      # if not (because of the algorithm or because of blacklists) we need to
      # test all possible subsets of the intersection of the markov blankets.
      if (!is.null(dsep.set) && !blacklisted.arc) {

        el = dsep.set[[which(sapply(dsep.set, function(x) setequal(x$arc, c(y, z))))]]

        # an unshielded triplet is a v-structure if and only if the central node
        # (x) is not part of the d-separating set of the other nodes (y and z).
        if (x %!in% el$dsep.set) {

          if (debug)
            cat("    @ detected v-structure", y, "->", x, "<-", z, "from d-separating set.\n")

          vs = rbind(vs, data.frame(max_a = el$p.value, y, x, z, stringsAsFactors = FALSE))

        }#THEN

      }#THEN
      else {

        # choose the smallest of mb(y) - {x,z} and mb(z) - {x,y} to cut down
        # the number of subsets to test.
        sx = smaller(setdiff(mb[[y]][['mb']], c(x, z)),
                     setdiff(mb[[z]][['mb']], c(x, y)))

        if (debug)
          cat("    > chosen d-separating set: '", sx, "'\n")

        a = allsubs.test(x = y, y = z, fixed = x, sx = sx, data = data,
              test = test, B = B, alpha = alpha, max = min(max.sx, length(sx)),
              complete = complete, debug = debug)

        if (a["p.value"] <= alpha) {

          if (debug)
            cat("    @ detected v-structure", y, "->", x, "<-", z, "\n")

          vs = rbind(vs, data.frame(max_a = a["max.p.value"], y, x, z, stringsAsFactors = FALSE))

        }#THEN

      }#ELSE

    }#FOR

    return(vs)

  }#VSTRUCT.CENTERED.ON

  sapply(nodes, vstruct.centered.on, mb = mb, data = data,
    dsep.set = attr(mb, "dsep.set"), simplify = FALSE)

}#VSTRUCT.DETECT

# include v-structures in the network, setting the corresponding arc directions.
vstruct.apply = function(arcs, vs, nodes, debug = FALSE) {

  if (debug)
    cat("----------------------------------------------------------------\n")

  for (i in seq(nrow(vs))) {

    x = vs[i, "x"]
    y = vs[i, "y"]
    z = vs[i, "z"]
    max_a = vs[i, "max_a"]

    # check whether the network already includes conflicting v-structures.
    if (!(is.listed(arcs, c(y, x)) && is.listed(arcs, c(z, x)))) {

      if (debug)
        cat("* not applying v-structure", y, "->", x, "<-", z, "(", max_a, ")\n")

      warning("vstructure ", y, " -> ", x, " <- ", z, " is not applicable, ",
        "because one or both arcs are oriented in the opposite direction.")

      next

    }#THEN

    # tentatively add the arcs that make up the v-structure.
    temp = set.arc.direction(y, x, arcs)
    temp = set.arc.direction(z, x, temp)

    # check whether the network is acyclic.
    if (!is.acyclic(temp, nodes, directed = TRUE)) {

      if (debug)
        cat("* not applying v-structure", y, "->", x, "<-", z, "(", max_a, ")\n")

      warning("vstructure ", y, " -> ", x, " <- ", z, " is not applicable, ",
        "because one or both arcs introduce cycles in the graph.")

      next

    }#THEN

    if (debug)
      cat("* applying v-structure", y, "->", x, "<-", z, "(", max_a, ")\n")

    # save the updated arc set.
    arcs = temp

  }#FOR

  return(arcs)

}#VSTRUCT.APPLY

# emergency measures for markov blanket and neighbourhood recovery.
bn.recovery = function(bn, nodes, filter = "AND", mb = FALSE, debug = FALSE) {

  .Call(call_bn_recovery,
        bn = bn,
        mb = mb,
        filter = match(filter, c("OR", "AND")),
        debug = debug)

}#BN.RECOVERY

# explore the structure of the network using its arc set.
cache.structure = function(nodes, arcs, amat = NULL, debug = FALSE) {

  # rebuild the adjacency matrix only if it's not available.
  if (is.null(amat))
    amat = arcs2amat(arcs, nodes)

  .Call(call_cache_structure,
        nodes = nodes,
        amat = amat,
        debug = debug)

}#CACHE.STRUCTURE

# explore the structure of the neighbourhood of a target node.
cache.partial.structure = function(nodes, target, arcs, amat = NULL,
    debug = FALSE) {

  # rebuild the adjacency matrix only if it's not available.
  if (is.null(amat))
    amat = arcs2amat(arcs, nodes)

  .Call(call_cache_partial_structure,
        nodes = nodes,
        target = target,
        amat = amat,
        debug = debug)

}#CACHE.PARTIAL.STRUCTURE
