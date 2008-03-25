
choose.direction.test = function(x, arc, data, test, alpha, debug) {

  nodes = names(x$nodes)
  amat = arcs2amat(x$arcs, nodes)

  # you can't help but notice nodes connected by undirected arcs are
  # included, too? wonder why?
  # because if they, too, are parents of the node to be tested
  # they _do_ belong there; if they are not, the node distribution
  # does not depend on them so they are largely irrelevant.

  parents1 = parents.backend(x$arcs, arc[2], undirected = TRUE)
  a1 = conditional.test(arc[1], arc[2],
        parents1[parents1 != arc[1]],
        data = data, test = test)

  parents2 = parents.backend(x$arcs, arc[1], undirected = TRUE)
  a2 = conditional.test(arc[2], arc[1],
        parents2[parents2 != arc[2]],
        data = data, test = test)

  if (debug) {

    cat("  > testing", arc[1], "->", arc[2], "with conditioning set '",
      parents1[parents1 != arc[1]], "'.\n")
    cat("    > p-value is", a1, ".\n")

    cat("  > testing", arc[2], "->", arc[1], "with conditioning set '",
      parents2[parents2 != arc[2]], "'.\n")
    cat("    > p-value is", a2, ".\n")

  }#THEN

  choose = function(a, b, x, arc) {

    cycles = has.path(arc[2], arc[1], nodes, amat, exclude.direct = TRUE)

    if (isTRUE(all.equal(as.numeric(a), as.numeric(b)))) {

      if (debug) cat("  @ nothing to do, same p-value.\n")
      return(x)

    }#THEN
    else if (a < b) {

      if (a < alpha) {

        if (cycles) {

          if (debug) {

            cat("  > adding", arc[1], "->", arc[2], "creates cycles!.\n")

          }#THEN

          # if one arc creates cycles, try the other one.
          if (b < alpha)
            choose(a = b, b = a, x = x, arc = arc[c(2, 1)])
          else {

            if (debug) cat("  > arc", arc[2], "->", arc[1], "isn't good, either.\n")
            return(x)

          }#ELSE

        }#THEN
        else {

          if (debug) cat("  @ arc", arc[1], "->", arc[2], "is better .\n")
          x$arcs = set.arc.direction(arc[1], arc[2], x$arcs)
          x$nodes = cache.structure(names(x$nodes), x$arcs)

          return(x)

        }#ELSE

      }#THEN
      else {

        # both tests are over the alpha threshold.
        if (debug) cat("  @ nothing to do, both p-values greater than", alpha, ".\n")
        return(x)

      }#ELSE

    }#THEN
    else if (b < a) {

      choose (a = b, b = a, x, arc = arc[c(2, 1)])

    }#THEN

  }#CHOOSE

  choose(a1, a2, x, arc)

}#CHOOSE.DIRECTION.TEST

choose.direction.score = function(x, data, arc, score, extra.args, debug) {

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
  if (score == "k2") {

    reference.score = sapply(arc, dirichlet.node, x = x, data = data)

  }#THEN
  else if (score %in% c("bde", "dir")) {

    reference.score = sapply(arc, dirichlet.node, x = x, data = data,
                        imaginary.sample.size = extra.args$iss)

  }#THEN
  else if (score %in% c("lik", "loglik")) {

    reference.score = sapply(arc, loglik.node, x = x, data = data)
    if (score == "lik") reference.score = exp(reference.score)

  }#THEN
  else if (score %in% c("aic", "bic")) {

    reference.score = sapply(arc, aic.node, x = x, data = data, k = extra.args$k)

  }#THEN
  else if (score == "bge") {

    reference.score = sapply(arc, bge.node, x = x, data = data,
                        imaginary.sample.size = extra.args$iss,
                        phi = extra.args$phi)

  }#THEN

  # compare the scores of the two networks.
  better1 = score.delta(arc = arc, network = x, data = data,
              score = score, score.delta = 0,
              reference.score = reference.score, op = "set",
              extra = extra.args, debug = debug)

  better2 = score.delta(arc = arc[c(2,1)], network = x, data = data,
              score = score, score.delta = 0,
              reference.score = reference.score, op = "set",
              extra = extra.args, debug = debug)

  if (debug) {

    cat("  > initial score for node", arc[1], "is", reference.score[1], ".\n")
    cat("  > initial score for node", arc[2], "is", reference.score[2], ".\n")
    cat("  > score delta for arc", arc[1], "->", arc[2], "is", better1$delta, ".\n")
    cat("  > score delta for arc", arc[2], "->", arc[1], "is", better2$delta, ".\n")

  }#THEN

  choose = function(a, b, x, arc) {

    cycles = has.path(arc[2], arc[1], nodes, amat)

    if (isTRUE(all.equal(as.numeric(better1$delta), as.numeric(better2$delta)))) {

      if (debug) cat("  @ nothing to do, same score delta.\n")
      return(x2)

    }#THEN
    else if (a$delta > b$delta) {

      if (a$bool) {

        if (cycles) {

          if (debug) {

            cat("  > adding", arc[1], "->", arc[2], "creates cycles!.\n")

          }#THEN

          # if one arc creates cycles, try the other one.
          if (b$bool)
            choose(a = b, b = a, x = x, arc = arc[c(2, 1)])
          else {

            if (debug) cat("  > arc", arc[2], "->", arc[1], "isn't good, either.\n")
            return(x2)

          }#ELSE

        }#THEN
        else {

          if (debug) cat("  @ arc", arc[1], "->", arc[2], "is better .\n")
          x$arcs = set.arc.direction(arc[1], arc[2], x$arcs)
          x$nodes = cache.structure(names(x$nodes), x$arcs)

          return(x)

        }#ELSE

      }#THEN
      else {

        # both tests are over the alpha threshold.
        if (debug) cat("  @ nothing to do, both score delta are negative.\n")
        return(x2)

      }#ELSE

    }#THEN
    else if (b$delta > a$delta) {

      choose (a = b, b = a, x, arc = arc[c(2, 1)])

    }#THEN

  }#CHOOSE

  choose(a = better1, b = better2, x = x, arc = arc)

}#CHOOSE.DIRECTION.SCORE
