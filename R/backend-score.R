
# compare two network scores in an efficient way.
score.delta = function(arc, network, data, score, score.delta,
    reference.score, op, extra, debug = FALSE) {

  # do a minimal update to the network structure.
  fake = .Call("score_delta_helper",
          net = network,
          arc = arc,
          operator = op,
          PACKAGE = "bnlearn")

  if (op == "reverse") {

    # compute the updated score contributions of the nodes involved.
    new.score = per.node.score(network = fake, score = score,
                        nodes = arc, extra.args = extra, data = data)

    # update the test counter.
    assign(".test.counter", get(".test.counter", envir = .GlobalEnv) + 2,
      envir = .GlobalEnv)

    # compare the network scores, minus numeric tolerance for better score
    # equivalence detection.
    if (isTRUE(all.equal(sum(new.score), sum(reference.score[arc]))))
      retval = 0
    else
      retval = sum(new.score) - sum(reference.score[arc])

  }#THEN
  else {

    # compute the updated score contributions of arc[2].
    new.score = per.node.score(network = fake, score = score,
                        nodes = arc[2], extra.args = extra, data = data)

    # update the test counter.
    assign(".test.counter", get(".test.counter", envir = .GlobalEnv) + 1,
      envir = .GlobalEnv)

    # compare the network scores.
    retval = new.score - reference.score[arc[2]]

  }#ELSE

  if (debug)
    cat("    > delta between scores for nodes", arc, "is", retval, ".\n")

  return(list(bool = (retval > score.delta + sqrt(.Machine$double.eps)),
    delta = retval, updates = new.score))

}#SCORE.DELTA

# create a data frame or an adjacency matrix containing the arcs to be added.
arcs.to.be.added = function(amat, nodes, blacklist = NULL, whitelist = NULL,
    arcs = TRUE) {

  .Call("hc_to_be_added",
        arcs = amat,
        blacklist = blacklist,
        whitelist = whitelist,
        nodes = nodes,
        convert = arcs,
        PACKAGE = "bnlearn")

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
arcs.to.be.reversed = function(arcs, blacklist) {

  if (!is.null(blacklist))
    return(arcs[!which.listed(arcs[, c(2, 1), drop = FALSE], blacklist), , drop = FALSE])
  else
    return(arcs)

}#ARCS.TO.BE.REVERSED

