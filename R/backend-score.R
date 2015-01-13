
# compare two network scores in an efficient way.
score.delta = function(arc, network, data, score, score.delta,
    reference.score, op, extra, decomposable = TRUE) {

  .Call("score_delta",
        arc = arc,
        network = network,
        data = data,
        score = score,
        score.delta = score.delta,
        reference.score = reference.score,
        op = op,
        extra = extra,
        decomposable = decomposable)

}#SCORE.DELTA

# create a data frame or an adjacency matrix containing the arcs to be added.
arcs.to.be.added = function(amat, nodes, blacklist = NULL, whitelist = NULL,
    nparents = NULL, maxp = Inf, arcs = TRUE) {

  .Call("hc_to_be_added",
        arcs = amat,
        blacklist = blacklist,
        whitelist = whitelist,
        nparents = nparents,
        maxp = maxp,
        nodes = nodes,
        convert = arcs)

}#ARCS.TO.BE.ADDED

# create a data frame containing the arcs to be dropped:
# arcs                   arcs already in the graph.
# !is.listed(whitelist)  exclude whitelisted arcs.
arcs.to.be.dropped = function(arcs, whitelist) {

  if (!is.null(whitelist))
    return(arcs[!which.listed(arcs, whitelist), , drop = FALSE])
  else
    return(arcs)

}#ARCS.TO.BE.DROPPED

# create a data frame containing the arcs to be reversed:
arcs.to.be.reversed = function(arcs, blacklist, nparents, maxp = Inf) {

  if (!is.null(blacklist))
    arcs = arcs[!which.listed(arcs[, c(2, 1), drop = FALSE], blacklist), , drop = FALSE]

  if (!missing(nparents))
    arcs = arcs[nparents[arcs[, 1]] < maxp, , drop = FALSE]

  return(arcs)

}#ARCS.TO.BE.REVERSED

