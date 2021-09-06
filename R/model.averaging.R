
# model averaging for bootstrapped network structures.
averaged.network.backend = function(strength, threshold) {

  nodes = attr(strength, "nodes")
  e = empty.graph(nodes)

  # arcs with a strength of one should always be selected, regardless of
  # the threshold.
  significant = (strength$strength > threshold) | (strength$strength == 1)

  # filter also the direction if present in the bn.strength object.
  if ("direction" %in% names(strength))
    significant = significant & (strength$direction >= 0.5)

  # nothing to see, move along.
  if (!any(significant))
    return(e)

  candidate.arcs = as.matrix(strength[significant, c("from", "to"), drop = FALSE])

  if (all(which.undirected(candidate.arcs))) {

    # update the arcs of the network, no cycles.
    e$arcs = candidate.arcs

  }#THEN
  else {

    # update the arcs of the network, minding cycles.
    e$arcs = .Call(call_smart_network_averaging,
                   arcs = candidate.arcs,
                   nodes = nodes,
                   weights = strength$strength[significant])

  }#ELSE

  # update the network structure.
  e$nodes = cache.structure(nodes, arcs = e$arcs)
  # add back illegal arcs, so that cpdag() works correctly.
  if ("illegal" %in% names(attributes(strength)))
    e$learning$illegal = attr(strength, "illegal")

  return(e)

}#AVERAGED.NETWORK.BACKEND

# compute the significance threshold for Friedman's confidence.
threshold = function(strength, method = "l1") {

  # do not blow up with graphs with only 1 node.
  if (nrow(strength) == 0)
    return(0)

  e = ecdf(strength$strength)
  u = knots(e)

  if (method == "l1") {

    norm = function(p)
      sum( diff(unique(c(0, u, 1))) * abs(e(unique(c(0, u[u < 1]))) - p))

  }#THEN

  p0 = optimize(f = norm, interval = c(0, 1))$minimum

  # double-check the boundaries, they are legal solutions but optimize() does
  # not check them.
  if (norm(1) < norm(p0))
    p0 = 1
  if (norm(0) < norm(p0))
    p0 = 0

  quantile(strength$strength, p0, type = 1, names = FALSE)

}#THRESHOLD

