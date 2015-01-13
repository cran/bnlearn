
is.legal.arc = function(arc, test, data) {

  if (test %in% c(available.mixedcg.tests, available.mixedcg.scores))
    if (is(data[, arc[1]], "numeric") && is(data[, arc[2]], "factor"))
      return(FALSE)

  return(TRUE)

}#IS.LEGAL.ARC

choose.direction.decide = function(x, arc, a, b, t, criterion, debug) {

  nodes = names(x$nodes)
  updated.nodes = character(0)

  # check whether any of the two directions causes cycles in the graph.
  amat = arcs2amat(x$arcs, nodes)
  cycles1 = has.path(arc[2], arc[1], nodes, amat, exclude.direct = TRUE)
  cycles2 = has.path(arc[1], arc[2], nodes, amat, exclude.direct = TRUE)

  if (criterion == "score") {

    good.enough = function(a) a > t
    better = function(a, b) a > b
    same = function(a, b) isTRUE(all.equal(as.numeric(a), as.numeric(b)))
    err.same = "  @ nothing to do, same score delta modulo precision tolerance.\n"
    err.neg = "  @ nothing to do, both score deltas are negative.\n"

  }#THEN
  else if (criterion == "boot") {

    good.enough = function(a) TRUE
    better = function(a, b) a$direction > b$direction
    same = function(a, b) isTRUE(all.equal(as.numeric(a$direction), as.numeric(b$direction)))
    err.same = "  @ nothing to do, same confidence modulo precision tolerance.\n"
    err.neg = paste0("  @ nothing to do, confidence is less than ", t, ".\n")

    if (a$strength < t) {

      if (debug)
         cat(err.neg)

      return(x)

    }#THEN

  }#THEN
  else if (criterion == "test") {

    good.enough = function(a) !is.na(a) && (a < t)
    better = function(a, b) ifelse(is.na(a), Inf, a) < ifelse(is.na(b), Inf, b)
    same = function(a, b) isTRUE(all.equal(as.numeric(a), as.numeric(b)))
    err.same = "  @ nothing to do, same p-value modulo precision tolerance.\n"
    err.neg = paste0("  @ nothing to do, p-value is less than ", t, ".\n")

  }#THEN

  if (cycles1 && cycles2) {

    # bothe A -> B and B -> A introduce cycles (should not be here).
    if (debug)
      cat("  @ nothing to do, both arc create cycles.\n")

  }#THEN
  else if (cycles1 && !cycles2) {

    # A -> B introduces cycles, B -> A does not.
    if (debug)
      cat("  > adding", arc[1], "->", arc[2], "creates cycles!.\n")

    if (good.enough(b)) {

      if (debug)
        cat("  @ arc", arc[2], "->", arc[1], "is better .\n")

       # update the arc set.
       x$arcs = set.arc.direction(arc[2], arc[1], x$arcs)
       # check which nodes have to be updated.
       updated.nodes = unique(c(arc, x$nodes[[arc[1]]]$mb, x$nodes[[arc[2]]]$mb))

    }#THEN
    else {

      if (debug)
        cat("  > arc", arc[2], "->", arc[1], "isn't good, either.\n")

    }#ELSE

  }#THEN
  else if (cycles2 && !cycles1) {

    # B -> A introduces cycles, A -> B does not.
    if (debug)
      cat("  > adding", arc[2], "->", arc[1], "creates cycles!.\n")

    if (good.enough(a)) {

      if (debug)
        cat("  @ arc", arc[1], "->", arc[2], "is better .\n")

       # update the arc set.
       x$arcs = set.arc.direction(arc[1], arc[2], x$arcs)
       # check which nodes have to be updated.
       updated.nodes = unique(c(arc, x$nodes[[arc[1]]]$mb, x$nodes[[arc[2]]]$mb))

    }#THEN
    else {

      if (debug)
        cat("  > arc", arc[1], "->", arc[2], "isn't good, either.\n")

    }#ELSE

  }#THEN
  else if (!cycles1 && !cycles2) {

    # neither A -> B nor B -> A introduce cycles.
    if (same(a, b)) {

      if (debug)
        cat(err.same)

    }#THEN
    else if (!good.enough(a) && !good.enough(b)) {

      if (debug)
        cat(err.neg)

    }#THEN
    else {

      id = if (better(a, b)) 1:2 else 2:1

      if (debug)
        cat("  @ arc", arc[id[1]], "->", arc[id[2]], "is better .\n")

       # update the arc set.
       x$arcs = set.arc.direction(arc[id[1]], arc[id[2]], x$arcs)
       # check which nodes have to be updated.
       updated.nodes = unique(c(arc, x$nodes[[arc[1]]]$mb, x$nodes[[arc[2]]]$mb))

    }#ELSE

  }#THEN

  # update the chosen nodes.
  for (node in updated.nodes)
    x$nodes[[node]] = cache.partial.structure(nodes, target = node,
      arcs = x$arcs, debug = FALSE)

  return(x)

}#CHOOSE.DIRECTION.DECIDE

choose.direction.test = function(x, arc, data, test, alpha, B, debug = FALSE) {

  # you can't help but notice nodes connected by undirected arcs are
  # included, too? wonder why?
  # because if they, too, are parents of the node to be tested
  # they _do_ belong there; if they are not, the node distribution
  # does not depend on them so they are largely irrelevant.
  choose.direction.test.pvalue = function(arc) {

    if (!is.legal.arc(arc, test, data)) {

      if (debug)
        cat("  > arc", arc[1], "->", arc[2], "is not a valid arc for this network.\n")

       return(NaN)

    }#THEN

    parents = parents.backend(x$arcs, arc[2], undirected = TRUE)
    a = indep.test(arc[1], arc[2], parents[parents != arc[1]], data = data,
          test = test, B = B, alpha = alpha)

    if (debug) {

      cat("  > testing", arc[1], "->", arc[2], "with conditioning set '",
        parents[parents != arc[1]], "'.\n")
      cat("    > p-value is", a, ".\n")

    }#THEN

    return(a)

  }#CHOOSE.DIRECTION.TEST.PVALUE

  a1 = choose.direction.test.pvalue(arc)
  a2 = choose.direction.test.pvalue(arc[2:1])

  choose.direction.decide(x = x, arc = arc, a = a1, b = a2, t = alpha,
    criterion = "test", debug = debug)

}#CHOOSE.DIRECTION.TEST

choose.direction.score = function(x, data, arc, score, extra.args, debug = FALSE) {

  # do a backup copy of the network structure.
  x2 = x

  # drop any existing arc between arc["from"] and arc["to"].
  x$arcs = drop.arc.backend(x$arcs, arc)

  # you can't help but notice nodes connected by undirected arcs are
  # included, too? wonder why?
  # because if they, too, are parents of the node to be tested
  # they _do_ belong there; if they are not, the node distribution
  # does not depend on them so they are largely irrelevant.

  x$nodes[[arc[1]]]$parents = parents.backend(x$arcs, arc[1], undirected = TRUE)
  x$nodes[[arc[2]]]$parents = parents.backend(x$arcs, arc[2], undirected = TRUE)
  x$nodes[[arc[1]]]$parents = x$nodes[[arc[1]]]$parents[x$nodes[[arc[1]]]$parents != arc[2]]
  x$nodes[[arc[2]]]$parents = x$nodes[[arc[2]]]$parents[x$nodes[[arc[2]]]$parents != arc[1]]

  # compute the initial score of the nodes involved.
  reference.score = per.node.score(network = x, score = score,
                      targets = arc, extra.args = extra.args, data = data)
  # check whether the score is decomposable.
  decomp = is.score.decomposable(score, names(x$nodes), extra.args)

  # compare the scores of the two networks.
  choose.direction.score.delta = function(arc) {

    better = score.delta(arc = arc, network = x, data = data,
               score = score, score.delta = 0,
               reference.score = reference.score, op = "set",
               extra = extra.args, decomposable = decomp)

    if (debug) {

      cat("  > initial score for node", arc[1], "is", reference.score[arc[1]], ".\n")
      if (is.legal.arc(arc, score, data))
        cat("  > score delta for arc", arc[1], "->", arc[2], "is", better$delta, ".\n")
      else
        cat("  > arc", arc[1], "->", arc[2], "is not a valid arc for this network.\n")

    }#THEN

    return(better)

  }#CHOOSE.DIRECTION.SCORE.DELTA

  a = choose.direction.score.delta(arc)
  b = choose.direction.score.delta(arc[2:1])

  choose.direction.decide(x = x2, arc = arc, a = a$delta, b = b$delta, t = 0,
    criterion = "score", debug = debug)

}#CHOOSE.DIRECTION.SCORE

choose.direction.boot = function(x, data, arc, extra.args, algorithm,
    algorithm.args, cpdag = TRUE, debug = FALSE) {

  # build a separate arc set with the two directions of the arc.
  m = matrix(c(arc, rev(arc)), ncol = 2, byrow = TRUE,
        dimnames = list(c(), c("from", "to")))

  # compute the respective bootstrap strength/direction
  res = arc.strength.boot(data = data, R = extra.args$R, m = extra.args$m,
          algorithm = algorithm, algorithm.args = algorithm.args,
          arcs = m, cpdag = cpdag, debug = FALSE)


  if (debug) {

    cat("  > testing", arc[1], "->", arc[2], "\n")
    cat("    > bootstrap probability of an arc between", arc[1], "and", arc[2],
      "is", res[1, "strength"], ".\n")
    if (is.legal.arc(arc, "loglik-cg", data))
      cat("    > direction confidence for arc", arc[1], "->", arc[2], "is",
        res[1, "direction"], ".\n")
    else
      cat("    > arc", arc[1], "->", arc[2], "is not a valid arc for this network.\n")
    if (is.legal.arc(arc[2:1], "loglik-cg", data))
      cat("    > direction confidence for arc", arc[2], "->", arc[1], "is",
        res[2, "direction"], ".\n")
    else
      cat("    > arc", arc[2], "->", arc[1], "is not a valid arc for this network.\n")

  }#THEN

  choose.direction.decide(x = x, arc = arc, a = res[1, ], b = res[2, ], t = 0.5,
    criterion = "boot", debug = debug)

}#CHOOSE.DIRECTION.BOOT

