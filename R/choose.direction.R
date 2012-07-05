
choose.direction.test = function(x, arc, data, test, alpha, B, debug = FALSE) {

  nodes = names(x$nodes)
  amat = arcs2amat(x$arcs, nodes)

  # you can't help but notice nodes connected by undirected arcs are
  # included, too? wonder why?
  # because if they, too, are parents of the node to be tested
  # they _do_ belong there; if they are not, the node distribution
  # does not depend on them so they are largely irrelevant.

  parents1 = parents.backend(x$arcs, arc[2], undirected = TRUE)
  a1 = conditional.test(arc[1], arc[2], parents1[parents1 != arc[1]],
        data = data, test = test, B = B, alpha = alpha)

  parents2 = parents.backend(x$arcs, arc[1], undirected = TRUE)
  a2 = conditional.test(arc[2], arc[1], parents2[parents2 != arc[2]],
        data = data, test = test, B = B, alpha = alpha)

  if (debug) {

    cat("  > testing", arc[1], "->", arc[2], "with conditioning set '",
      parents1[parents1 != arc[1]], "'.\n")
    cat("    > p-value is", a1, ".\n")

    cat("  > testing", arc[2], "->", arc[1], "with conditioning set '",
      parents2[parents2 != arc[2]], "'.\n")
    cat("    > p-value is", a2, ".\n")

  }#THEN

  choose = function(a, b, x, arc, recurse = FALSE) {

    cycles = has.path(arc[2], arc[1], nodes, amat, exclude.direct = TRUE)

    if (isTRUE(all.equal(as.numeric(a), as.numeric(b)))) {

      if (debug)
        cat("  @ nothing to do, same p-value.\n")

      return(x)

    }#THEN
    else if (a < b) {

      if (a < alpha) {

        if (cycles) {

          if (debug)
            cat("  > adding", arc[1], "->", arc[2], "creates cycles!.\n")

          # if one arc creates cycles, try the other one.
          if ((b < alpha) && recurse)
            choose(a = b, b = a, x = x, arc = arc[c(2, 1)], recurse = FALSE)
          else {

            if (debug)
              cat("  > arc", arc[2], "->", arc[1], "isn't good, either.\n")

            return(x)

          }#ELSE

        }#THEN
        else {

          if (debug)
            cat("  @ arc", arc[1], "->", arc[2], "is better.\n")

          # update the arc set.
          x$arcs = set.arc.direction(arc[1], arc[2], x$arcs)
          # check which nodes have to be updated.
          updated.nodes = unique(c(arc, x$nodes[[arc[1]]]$mb, x$nodes[[arc[2]]]$mb))
          # update the chosen nodes.
          for (node in updated.nodes)
            x$nodes[[node]] = cache.partial.structure(names(x$nodes),
              target = node, arcs = x$arcs, debug = FALSE)

          return(x)

        }#ELSE

      }#THEN
      else {

        # both tests are over the alpha threshold.
        if (debug)
          cat("  @ nothing to do, both p-values greater than", alpha, ".\n")

        return(x)

      }#ELSE

    }#THEN
    else if (b < a) {

      choose (a = b, b = a, x, arc = arc[c(2, 1)])

    }#THEN

  }#CHOOSE

  choose(a1, a2, x, arc, recurse = TRUE)

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

  # cache node names and the adjacency matrix.
  nodes = names(x$nodes)
  amat = arcs2amat(x$arcs, nodes)

  # compute the initial score of the nodes involved.
  reference.score = per.node.score(network = x, score = score,
                      nodes = arc, extra.args = extra.args, data = data)

  # compare the scores of the two networks.
  better1 = score.delta(arc = arc, network = x, data = data,
              score = score, score.delta = 0,
              reference.score = reference.score, op = "set",
              extra = extra.args)

  better2 = score.delta(arc = arc[c(2, 1)], network = x, data = data,
              score = score, score.delta = 0,
              reference.score = reference.score, op = "set",
              extra = extra.args)

  if (debug) {

    cat("  > initial score for node", arc[1], "is", reference.score[1], ".\n")
    cat("  > initial score for node", arc[2], "is", reference.score[2], ".\n")
    cat("  > score delta for arc", arc[1], "->", arc[2], "is", better1$delta, ".\n")
    cat("  > score delta for arc", arc[2], "->", arc[1], "is", better2$delta, ".\n")

  }#THEN

  choose = function(a, b, x, arc, recurse = FALSE) {

    cycles = has.path(arc[2], arc[1], nodes, amat, exclude.direct = TRUE)

    if (isTRUE(all.equal(as.numeric(better1$delta), as.numeric(better2$delta)))) {

      if (debug)
        cat("  @ nothing to do, same score delta.\n")

      return(x2)

    }#THEN
    else if (a$delta > b$delta) {

      if (a$bool) {

        if (cycles) {

          if (debug)
            cat("  > adding", arc[1], "->", arc[2], "creates cycles!.\n")

          # if one arc creates cycles, try the other one.
          if ((b$bool) & recurse)
            choose(a = b, b = a, x = x, arc = arc[c(2, 1)], recurse = FALSE)
          else {

            if (debug)
              cat("  > arc", arc[2], "->", arc[1], "isn't good, either.\n")

            return(x2)

          }#ELSE

        }#THEN
        else {

          if (debug)
            cat("  @ arc", arc[1], "->", arc[2], "is better .\n")

          # update the arc set.
          x$arcs = set.arc.direction(arc[1], arc[2], x$arcs)
          # check which nodes have to be updated.
          updated.nodes = unique(c(arc, x$nodes[[arc[1]]]$mb, x$nodes[[arc[2]]]$mb))
          # update the chosen nodes.
          for (node in updated.nodes)
            x$nodes[[node]] = cache.partial.structure(names(x$nodes),
              target = node, arcs = x$arcs, debug = FALSE)

          return(x)

        }#ELSE

      }#THEN
      else {

        # both tests are over the alpha threshold.
        if (debug)
          cat("  @ nothing to do, both score delta are negative.\n")

        return(x2)

      }#ELSE

    }#THEN
    else if (b$delta > a$delta) {

      choose (a = b, b = a, x, arc = arc[c(2, 1)])

    }#THEN

  }#CHOOSE

  choose(a = better1, b = better2, x = x, arc = arc, recurse = TRUE)

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
     cat("    > bootstrap probability of an arc between", arc[1], "and", arc[2], "is", res[1, "strength"], ".\n")
     cat("    > direction confidence for arc", arc[1], "->", arc[2], "is", res[1, "direction"], ".\n")
     cat("    > direction confidence for arc", arc[2], "->", arc[1], "is", res[2, "direction"], ".\n")

  }#THEN

  if (res[1, "strength"] < 0.5) {

    if (debug)
       cat("  @ nothing to do, bootstrap probability is less than 0.50.\n")

    return(x)

  }#THEN

  # cache node names and the adjacency matrix.
  nodes = names(x$nodes)
  updated.nodes = character(0)
  amat = arcs2amat(x$arcs, nodes)

  # check whether any of the two directions causes cycles in the graph.
  cycles1 = has.path(arc[2], arc[1], nodes, amat, exclude.direct = TRUE)
  cycles2 = has.path(arc[1], arc[2], nodes, amat, exclude.direct = TRUE)

  if (cycles1 && cycles2) {

    # bothe A -> B and B -> A introduce cycles (should not be here).
    if (debug)
      cat("  @ nothing to do, both arc create cycles.\n")

  }#THEN
  else if (cycles1 && !cycles2) {

    # A -> B introduces cycles, B -> A does not.
    if (debug)
      cat("  > adding", arc[1], "->", arc[2], "creates cycles!.\n")

    if (res[2, "direction"] > res[1, "direction"]) {

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

    if (res[1, "direction"] > res[2, "direction"]) {

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
    if (res[1, "direction"] > res[2, "direction"]) {

      if (debug)
        cat("  @ arc", arc[1], "->", arc[2], "is better .\n")

       # update the arc set.
       x$arcs = set.arc.direction(arc[1], arc[2], x$arcs)
       # check which nodes have to be updated.
       updated.nodes = unique(c(arc, x$nodes[[arc[1]]]$mb, x$nodes[[arc[2]]]$mb))

    }#THEN
    else {

      if (debug)
        cat("  @ arc", arc[2], "->", arc[1], "is better .\n")

       # update the arc set.
       x$arcs = set.arc.direction(arc[2], arc[1], x$arcs)
       # check which nodes have to be updated.
       updated.nodes = unique(c(arc, x$nodes[[arc[1]]]$mb, x$nodes[[arc[2]]]$mb))

    }#ELSE

  }#THEN

  # update the chosen nodes.
  for (node in updated.nodes)
    x$nodes[[node]] = cache.partial.structure(nodes, target = node,
      arcs = x$arcs, debug = FALSE)

  return(x)

}#CHOOSE.DIRECTION.BOOT

