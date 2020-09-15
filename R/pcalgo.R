
pc.stable.backend = function(x, cluster = NULL, whitelist, blacklist, test,
  alpha, B, max.sx = ncol(x), complete, debug = FALSE) {

  nodes = names(x)
  nnodes = length(nodes)
  mb = structure(vector(length(nodes), mode = "list"), names = nodes)
  skeleton = subsets(nodes, 2)
  node.pairs =
    apply(skeleton, 1, function(x) list(arc = x, max.adjacent = nnodes - 1))
  nbr.size = rep(nnodes - 1, length(node.pairs))

  # set the size of the largest conditioning set using either the given limit
  # or the number of variables in the data, whichever is lower.
  max.dsep.size = min(max.sx, length(nodes) - 2)

  if (debug && max.dsep.size < length(nodes) - 2)
    cat("@ limiting conditioning sets to", max.dsep.size, "nodes.\n")

  # find out which nodes are adjacent.
  for (dsep.size in seq(from = 0, to = max.dsep.size)) {

    # perform the conditional independence tests.
    node.pairs[dsep.size <= nbr.size] =
      smartSapply(cluster, node.pairs[dsep.size <= nbr.size], pc.heuristic,
        data = x, alpha = alpha, B = B, whitelist = whitelist,
        blacklist = blacklist, test = test, skeleton = skeleton,
        dsep.size = dsep.size, complete, debug = debug)

    # find out which undirected arcs are still present.
    arcs.still.present = lapply(node.pairs, function(x) {

      if (x$p.value < alpha)
        return(x$arc)
      else
        return(NULL)

    })

    # update the skeleton.
    skeleton = do.call(rbind, arcs.still.present)

    # count how many nodes (at most) are adjacent to each of the endpoints.
    nbr.size = sapply(node.pairs, `[[`, "max.adjacent")

    # if that number is smaller than the current size of the d-separation set,
    # there are no valid conditioning sets to test.
    if (all(nbr.size <= dsep.size))
      break

    if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("* remaining arcs:\n")
      print(arcs.rbind(skeleton, skeleton, reverse2 = TRUE))

    }#THEN

  }#FOR

  # start encoding the bn object.
  skeleton = cache.structure(nodes, arcs.rbind(skeleton, skeleton, reverse2 = TRUE))
  # attach the d-separating sets.
  attr(skeleton, "dsep.set") = node.pairs

  return(skeleton)

}#PC.STABLE.BACKEND

pc.heuristic = function(pair, data, alpha, B, whitelist, blacklist, test,
    skeleton, dsep.size, complete, debug = FALSE) {

  arc = pair$arc

  # check whether the arc is blacklisted in both directions (so that we do not
  # include it) or whitelisted in at least one direction (so that we include it).
  if (is.whitelisted(whitelist, arc, either = TRUE))
    return(list(arc = arc, p.value = 0, dsep.set = NULL, max.adjacent = 0))
  else if (is.blacklisted(blacklist, arc, both = TRUE))
    return(list(arc = arc, p.value = 1, dsep.set = NULL, max.adjacent = 0))

  # check whether the nodes are already d-separated.
  if (!is.null(pair$dsep.set))
    return(pair)

  # only nodes that are adjacent to the arc endpoints are investigated in the
  # search for a d-separating set.
  nbr1 = union(skeleton[skeleton[, 2] == arc[1], 1], skeleton[skeleton[, 1] == arc[1], 2])
  nbr1 = setdiff(nbr1, arc[2])
  nbr2 = union(skeleton[skeleton[, 2] == arc[2], 1], skeleton[skeleton[, 1] == arc[2], 2])
  nbr2 = setdiff(nbr2, arc[1])

  # not enough nodes to form a d-separating set of the given size.
  if ((length(nbr1) < dsep.size) && (length(nbr2) < dsep.size))
    return(list(arc = arc, p.value = pair$p.value, dsep.set = pair$dsep.set,
      max.adjacent = 0))

  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* investigating", arc[1], "-", arc[2] ,
        ", d-separating sets of size", dsep.size, ".\n")
    cat("  > neighbours of", arc[1], ":", nbr1, "\n")

  }#THEN

  if (length(nbr1) >= dsep.size) {

    a1 = allsubs.test(x = arc[1], y = arc[2], sx = nbr1, min = dsep.size,
           max = dsep.size, data = data, test = test, alpha = alpha, B = B,
           complete = complete, debug = debug)

    if (a1["p.value"] > alpha)
      return(list(arc = arc, p.value = a1["p.value"],
        dsep.set = attr(a1, "dsep.set"), max.adjacent = 0))

  }#THEN

  if (debug)
    cat("  > neighbours of", arc[2], ":", nbr2, "\n")

  # there are cases in which checking the neighbours of the second endpoint
  # is redundant:
  #   * if d-separating set is the empty set, because marginal tests are
  #       symmetric;
  #   * if the d-separating sets are of size 1 (just single nodes), then it is
  #       trivial not test them again when the same node is adjacent to both
  #       endpoints;
  #   * if nbr1 and nbr2 are identical, and thus produce the same set of tests.
  if (dsep.size == 1)
    nbr2 = setdiff(nbr2, nbr1)

  if ((length(nbr2) >= dsep.size) && (dsep.size > 0) && !setequal(nbr1, nbr2)) {

    a2 = allsubs.test(x = arc[2], y = arc[1], sx = nbr2, min = dsep.size,
           max = dsep.size, data = data, test = test, alpha = alpha, B = B,
           complete = complete, debug = debug)

    if (a2["p.value"] > alpha)
      return(list(arc = arc, p.value = a2["p.value"],
        dsep.set = attr(a2, "dsep.set"), max.adjacent = 0))

  }#THEN

  return(list(arc = arc, p.value = 0, dsep.set = NULL,
    max.adjacent = max(length(nbr1), length(nbr2))))

}#PC.HEURISTIC

