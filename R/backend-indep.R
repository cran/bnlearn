
# second prinple of CI algorithms: infer arc orientation from graph structure.
second.principle = function(x, cluster = NULL, mb, whitelist, blacklist,
  test, alpha, data, strict, direction, debug) {

  nodes = names(x)

  # build a list of the undirected arcs in the graph, using the neighbourhoods
  # detected in markov.blanket().
  arcs = mb2arcs(mb, nodes)

  # apply blacklist to the arc set.
  to.drop = !apply(arcs, 1, function(x){ is.blacklisted(blacklist, x) })
  arcs = arcs[to.drop, , drop = FALSE]

  # 3. [Orient Edges]
  # 3.1 detect v-structures.
  vs = do.call("rbind",
         vstruct.detect(nodes = nodes, arcs = arcs, mb = mb, data = x,
           alpha = alpha, test = test, debug = debug))
  rownames(vs) = NULL

  if (!is.null(vs)) {

    # 3.2 sort them in p-value order.
    vs = vs[order(vs[,"max_a"], decreasing = FALSE),]
    # 3.3 apply.
    arcs = vstruct.apply(arcs = arcs, vs = vs, strict = strict, debug = debug)

  }#THEN

  # 4. [Remove Cycles] and 5. [Reverse Edges]
  # thou shalt not create cycles in the graph, it's acyclic!
  arcs = orient.edges(arcs = arcs, nodes = nodes,
           whitelist = whitelist, blacklist = blacklist,
           cluster = cluster, debug = debug)

  # 6. [Propagate Directions]
  arcs = propagate.directions(arcs = arcs, nodes = nodes, debug = debug)

  # arcs whitelisted in one direction (i.e. a -> b but not b -> a) are
  # never checked for cycles, do it now.
  if (!is.acyclic(arcs = arcs, nodes = nodes))
    stop("the graph contains cycles, possibly because of whitelisted nodes.")

  # save the status of the learning algorithm.
  learning = list(whitelist = whitelist, blacklist = blacklist, 
    test = test, args = list(alpha = alpha),
    ntests = get(".test.counter", envir = .GlobalEnv))

  # EXTRA [ESP]
  if (direction)
    arcs = set.directions(arcs = arcs, data = x, test = test,
             alpha = alpha, cluster = cluster, debug = debug)

  list(learning = learning, nodes = cache.structure(nodes, arcs = arcs),
    arcs = arcs)

}#SECOND.PRINCIPLE

# propagate directions to the undirected arcs.
propagate.directions = function(arcs, nodes, debug) {

  # build the adjacency matrix.
  amat = arcs2amat(arcs, nodes)
  # ignore directed arcs.
  undirected.arcs = arcs[which.undirected(arcs),]

  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* propagating directions for the following undirected arcs:\n")

  }#THEN

  if (nrow(undirected.arcs) == 0) {

    if (debug) cat("  > nothing to do.\n")

  }#THEN
  else apply(undirected.arcs, 1, function(arc) {

    # if there is a path from arc[1] to arc[2] besides arc[1] -> arc[2] ...
    if (has.path(arc[1], arc[2], nodes, amat, exclude.direct = TRUE)) {

      # ... set the direction to arc[1] -> arc[2].
      assign("arcs", set.arc.direction(arc[1], arc[2], arcs),
        envir = sys.frame(-2))

      if (debug) {

        cat("  > there's a path from", arc[1], "to", arc[2], ".\n")
        cat("    > removing", arc[2], "->", arc[1], ".\n")

      }#THEN

    }#THEN
    else {

      if (debug)
        cat("  > no path from", arc[1], "to", arc[2], ".\n")

    }#ELSE

  })

  arcs

}#PROPAGATE.DIRECTIONS

# build the neighbourhood of a node from the markov blanket.
neighbour = function(x, mb, data, alpha, whitelist, blacklist,
  backtracking = NULL, test, debug) {

  # save a prisitine copy of the markov blanket.
  nbrhood = mb[[x]]

  # if the markov blanket is empty there's nothing to do.
  if (length(nbrhood) == 0) {

    if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("* markov blanket of", x, "is empty; the neighbourhood too.\n")

    }#THEN

    return(list(mb = nbrhood, nbr = nbrhood))

  }#THEN

  known.good = known.bad = c()
  blacklisted = nbrhood[sapply(nbrhood,
          function(y) { is.blacklisted(blacklist, c(x,y), both = TRUE) })]
  whitelisted = nbrhood[sapply(nbrhood,
          function(y) { is.whitelisted(whitelist, c(x,y), either = TRUE) })]

  # whitelisted nodes are included (arc orientation is irrelevant),
  # and blacklisted nodes are removed if both directed arcs are banned
  # and both are not in the whitelist.
  nbrhood = nbrhood[!(nbrhood %in% blacklisted)]
  nbrhood = unique(c(nbrhood, whitelisted))

  # use backtracking for a further screening of the nodes to be checked.
  if (!is.null(backtracking)) {

    # neighbourhood is symmetric, so X \in N(Y) <=> Y \in N(X)
    known.good = names(backtracking[backtracking])

    # and vice versa X \not\in N(Y) <=> Y \not\in N(X)
    known.bad = names(backtracking[!backtracking])

    # known.bad nodes are not to be checked for inclusion and/or used in
    # the subsets.
    nbrhood = nbrhood[!(nbrhood %in% known.bad)]

  }#THEN

  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* detecting neighbourhood of", x, ".\n")
    cat("  * blacklisted nodes: '", blacklisted, "'\n")
    cat("  * whitelisted nodes: '", whitelisted, "'\n")
    cat("  * starting with neighbourhood: '", nbrhood, "'\n")

    if (!is.null(backtracking)) {

      cat("  * known good (backtracking): '", known.good, "'.\n")
      cat("  * known bad (backtracking): '", known.bad, "'.\n")

    }#THEN

  }#THEN

  if (length(nbrhood) > 1) {

    nbr = function(y, x, mb, test) {

      k = 0
      a = 0

      if (debug)
        cat("  * checking node", y, "for neighbourhood.\n")

      # choose the smaller set of possible d-separating nodes.
      dsep.set = smaller(mb[[x]][mb[[x]] != y], mb[[y]][mb[[y]] != x])
      if (debug)
        cat("    > dsep.set = '", dsep.set, "'\n")

      repeat {

        # create all the possible subsets of the markov blanket (excluding
        # the node to be tested for exclusion) of size k.
        dsep.subsets = subsets(length(dsep.set), k, dsep.set)

        for (s in 1:nrow(dsep.subsets)) {

          if (debug)
            cat("    > trying conditioning subset '", dsep.subsets[s,], "'.\n")

          a = conditional.test(x, y, dsep.subsets[s,], data = data, test = test)
          if (a > alpha) {

            if (debug)
              cat("    > node", y, "is not a neighbour of", x, ". ( p-value:", a, ")\n")

            # update the neighbourhood.
            assign("nbrhood", nbrhood[nbrhood != y], envir = sys.frame(-3) )

            break

          }#THEN
          else if (debug) {

              cat("    > node", y, "is still a neighbour of", x, ". ( p-value:", a, ")\n")

          }#THEN

        }#FOR

        # if the node was not removed from the markov blanket, increase
        # the size of the conditioning set and try again (if the size of
        # the markov blanket itself allows it).
        if ((a <= alpha) && (k < length(dsep.set))) {

          k = k+1

        }#THEN
        else {

          break

        }#ELSE

      }#REPEAT

    }#NBR

    # do not even try to remove whitelisted and backtracked (good) nodes.
    sapply(nbrhood[!(nbrhood %in% unique(c(whitelisted, known.good)))],
             nbr, x = x, mb =mb, test = test)

  }#THEN

  list (mb = mb[[x]], nbr = nbrhood)

}#NEIGHBOUR

# detect v-structures in the graph.
vstruct.detect = function(nodes, arcs, mb, data, alpha, test, debug) {

  vstruct.centered.on = function(x, mb, data) {

    if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("* v-structures centered on", x, ".\n")

    }#THEN

    tos = parents.backend(arcs, x, TRUE)

    if (length(tos) < 2) return(NULL)

    # build a list of possibile parents for the node x, i.e. all the subsets
    # of size 2 of the nodes connected to x by incoming arcs.
    tos.combs = subsets(length(tos), 2, tos)
    vs = NULL

    for (j in 1:nrow(tos.combs)) {

      y = tos.combs[j,1]
      z = tos.combs[j,2]

      if (debug)
        cat("  * checking", y, "->", x, "<-", z, "\n")

      # check there's no arc from y to z and vice versa.
      if(!is.listed(arcs, c(y, z)) &&
         !is.listed(arcs, c(z, y))) {

        mby = mb[[y]][['mb']]
        mbz = mb[[z]][['mb']]

        # compute mb(y) - {x,z} and mb(z) - {x,y}
        mby = mby[!(mby %in% c(x, z))]
        mbz = mbz[!(mbz %in% c(x, y))]

        # choose the smallest one to cut down the number of subsets to test.
        dsep.set = smaller(mby, mbz)

        if (debug)
          cat("    > chosen d-separating set: '", dsep.set, "'\n")

        k = 0
        max_a = a = 0

        repeat {

          dsep.subsets = subsets(length(dsep.set), k, dsep.set)

          for (s in 1:nrow(dsep.subsets)) {

            a = conditional.test(y, z, c(dsep.subsets[s,], x), data = data, test = test)
            if (debug)
              cat("    > testing", y, "vs", z, "given", c(dsep.subsets[s,], x), "(", a, ")\n")
            max_a = max(a, max_a)
            if ((a > alpha) && debug) {

              cat("    >", y, "and", z, "are independent given '", c(dsep.subsets[s,], x), "' (", a, ")\n")
              break

            }#THEN

          }#FOR

          if (a <= alpha) {

            if (k < length(dsep.set)) {

              k = k + 1

            }#THEN
            else {

              if (debug)
                cat("    @ detected v-structure", y, "->", x, "<-", z, "\n")
              vs = rbind(vs, data.frame(max_a, y, x, z, stringsAsFactors = FALSE))
              break

            }#ELSE

          }#THEN
          else {

            break

          }#ELSE

        }#REPEAT

      }#THEN

    }#FOR

    vs

  }#VSTRUCT.CENTERED.ON

  sapply(nodes, vstruct.centered.on, mb = mb, data = data)

}#VSTRUCT.DETECT

# apply v-structures to a graph.
vstruct.apply = function(arcs, vs, strict, debug) {

  if (debug)
    cat("----------------------------------------------------------------\n")

  apply(vs, 1, function(v) {

    # check the arcs of the v-strustures are still there.
    if(is.listed(arcs, v[c("y", "x")]) &&
       is.listed(arcs, v[c("z", "x")])) {

      if (debug)
        cat("* applying v-structure", v["y"], "->", v["x"], "<-", v["z"],
              "(", v["max_a"], ")\n")

      assign("arcs", set.arc.direction(v["y"], v["x"], arcs),
        envir = sys.frame(-2))
      assign("arcs", set.arc.direction(v["z"], v["x"], arcs),
        envir = sys.frame(-2))

    }#THEN
    else {

      if (debug) {

        cat("* not applying v-structure", v["y"], "->", v["x"], "<-", v["z"],
              "(", v["max_a"], ")\n")
        cat("  > adding arc", v["y"], "->", v["z"], "\n")
        cat("  > adding arc", v["z"], "->", v["y"], "\n")

      }#THEN

      if (strict)
        stop(paste("vstructure", v["y"], "->", v["x"], "<-", v["z"],
          "is not applicable, because one or both arcs are oriented",
          "in the opposite direction."))
      else
        warning(paste("vstructure", v["y"], "->", v["x"], "<-", v["z"],
          "is not applicable, because one or both arcs are oriented",
          "in the opposite direction."))

    }#ELSE

  })

  arcs

}#VSTRUCT.APPLY

# remove arcs from the graph to make it acyclic.
orient.edges = function(arcs, nodes, whitelist, blacklist, cluster, debug) {

  to.be.reversed = matrix(character(0), ncol = 2,
                     dimnames = list(c(), c("from", "to")))
  narcs = nrow(arcs)
  n = 0

  # remove undirected arcs from the whitelist; their direction is
  # left to the algorithm to choose.
  if (!is.null(whitelist) && (nrow(whitelist) > 0))
    whitelist = whitelist[!which.undirected(whitelist),, drop = FALSE]

  repeat {

    # cannot remove or reverse more arcs than you have.
    if (n > narcs)
      stop("too many iterations, probably would have gone on forever.")

    # check arcs between nodes which are both part of at least one cycle.
    in.cycles = is.acyclic.backend(arcs, nodes, directed = FALSE, return.nodes = TRUE)
    to.be.checked = (arcs[, "from"] %in% in.cycles) & (arcs[, "to"] %in% in.cycles)

    if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("* detecting cycles ...\n")
      cat("  > ignored nodes:", length(nodes) - length(in.cycles), "/", length(nodes), "\n")
      cat("  > ignored arcs:", length(which(!to.be.checked)), "/", nrow(arcs), "\n")
      print(arcs[!to.be.checked, , drop = FALSE])
      cat("  > checked arcs:",length(which(to.be.checked)), "/", nrow(arcs), "\n")

    }#THEN

    # compute the cycle counter of each (relevant) arc.
    cycles = cycle.counter(arcs = arcs[to.be.checked,, drop = FALSE],
              nodes = in.cycles, cluster = cluster)

    # do not check arcs whitelisted directed arcs.
    cycles = cycles[!which.listed(cycles[, 1:2], whitelist),, drop = FALSE]

    if (debug && (nrow(cycles) > 0)) print(cycles)

    # if there are no more arcs to check, break.
    if (nrow(cycles) == 0) break

    # remove the arc belonging to the highest number of cycles.
    arcs = arcs[!is.row.equal(arcs, cycles[1, c("from", "to")]),, drop = FALSE]

    # if the arc was already directed and its reverse is not blacklisted,
    # store it for future reversal.
    if (is.directed(cycles[1, c("from", "to")], arcs) &&
        !is.blacklisted(blacklist, cycles[1, c("to", "from")])) {

      to.be.reversed = rbind(to.be.reversed,
                         as.matrix(cycles[1, c("from", "to")]))

    }#THEN

    if (debug) {

      cat("  > arc", cycles[1, "from"], "->", cycles[1, "to"],
        "will be removed (", cycles[1, "cycles"], " cycles ).\n")

      cat("  > arcs scheduled for reversal are:\n")
      print(to.be.reversed)

    }#THEN

    n = n + 1

  }#REPEAT

  # add removed arcs, reversed and return.
  rbind(arcs, to.be.reversed[, c(2,1)])

}#ORIENT.EDGES

# count how many cycles each arc is part of.
cycle.counter = function(arcs, nodes, cluster, debug = FALSE) {

  # build the adjacency matrix.
  amat = arcs2amat(arcs, nodes)
  # cycle counter
  counter = data.frame(arcs, cycles = rep(0, nrow(arcs)),
    stringsAsFactors = FALSE)

  # if there are no arcs, there is nothing to do here.
  if (nrow(arcs) == 0)
    return(counter)

  if (!is.null(cluster)) {

    counter[, "cycles"] = parApply(cluster, arcs, 1, how.many.cycles,
      nodes = nodes, amat = amat, debug = FALSE)

  }#THEN
  else {

    counter[, "cycles"] = apply(arcs, 1, how.many.cycles, nodes = nodes,
      amat = amat, debug = debug)

  }#ELSE

  # return the sorted cycle counter.
  counter[order(counter[,"cycles"], decreasing = TRUE),, drop = FALSE]

}#CYCLE.COUNTER

# test undirected arcs in both direction to infer their orientation.
set.directions = function(arcs, data, test, alpha, cluster, debug) {

  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* infer arc direction.\n")

  }#THEN

  # get the undirected arcs.
  undirected.arcs = arcs[which.undirected(arcs), , drop = FALSE]

  if (nrow(undirected.arcs) == 0)
    return(arcs)

  if (is.null(cluster)) {

    a = apply(undirected.arcs, 1, function(arc) {

      # you can't help but notice nodes connecetd by undirected arcs are
      # included, too? wonder why?
      # because if they, too, are parents of the node to be tested
      # they _do_ belong there; if they are not, the node distribution
      # does not depend on them so they are largely irrelevant.
      parents = parents.backend(arcs, arc[2], TRUE)
      conditional.test(arc[1], arc[2],
          parents[parents != arc[1]],
          data = data, test = test)

    })

  }#THEN
  else {

     a = parApply(cluster, undirected.arcs, 1, function(arc, data, test) {

      # you can't help but notice nodes connecetd by undirected arcs are
      # included, too? wonder why?
      # because if they, too, are parents of the node to be tested
      # they _do_ belong there; if they are not, the node distribution
      # does not depend on them so they are largely irrelevant.
      parents = parents.backend(arcs, arc[2], TRUE)
      conditional.test(arc[1], arc[2],
          parents[parents != arc[1]],
          data = data, test = test)

    }, data = data, test = test)

  }#ELSE

  # build the lookup table.
  lookup.table = cbind(undirected.arcs, a)

  # sort arcs by increasing p-value.
  lookup.table = lookup.table[order(as.numeric(a)), ]

  repeat {

    # no more arcs, breaking the loop.
    if(nrow(lookup.table) == 0) break

    # get the most likely arc.
    da.best = lookup.table[1, ]

    if (debug) {

      cat("  > the results of the conditional independence tests are:\n")
      print(lookup.table)

    }#THEN

    # if both arcs have the same p-value, do not remove them from the arc set.
    if (da.best[3] != lookup.table[is.row.equal(lookup.table[, 1:2], da.best[c(2, 1)]), 3]) {

      # remove its reverse from the arc set.
      arcs = set.arc.direction(da.best[1], da.best[2], arcs)

      if (debug) {

        cat("  > ", da.best[1], "->", da.best[2], "has the lowest alpha (", da.best[3], ").\n")
        cat("  > ", da.best[2], "->", da.best[1], "will be removed.\n")

      }#THEN

    }#THEN
    else {

      if (debug) cat("  > it's a draw, no arc will be removed.\n")

    }#ELSE

    # remove both of them from the lookup table
    lookup.table = lookup.table[!is.row.equal(lookup.table[, 1:2], da.best[1:2]) &
                                !is.row.equal(lookup.table[, 1:2], da.best[c(2, 1)]), ]

  }#REPEAT

  arcs

}#SET.DIRECTIONS

# emergency measures for markov blanket and neighbourhood recovery.
bn.recovery = function(bn, nodes, strict, mb = FALSE, debug) {

  .Call("bn_recovery",
        bn = bn,
        strict = strict,
        mb = mb,
        debug = debug,
        PACKAGE = "bnlearn")

}#BN.RECOVERY

# explore the structure of the network using its arc set.
cache.structure = function(nodes, arcs, amat = NULL, debug = FALSE) {

  # rebuild the adjacency matrix only if it's not available
  if (is.null(amat))
    amat = arcs2amat(arcs, nodes)

  .Call("cache_structure",
        nodes = nodes,
        amat = amat,
        debug = debug,
        PACKAGE = "bnlearn")

}#CACHE.STRUCTURE

# explore the structure of the neighbourhood of a target node.
cache.partial.structure = function(nodes, target, arcs, amat = NULL,
    debug = FALSE) {

  # rebuild the adjacency matrix only if it's not available
  if (is.null(amat))
    amat = arcs2amat(arcs, nodes)

  .Call("cache_partial_structure",
        nodes = nodes,
        target = target,
        amat = amat,
        debug = debug,
        PACKAGE = "bnlearn")

}#CACHE.PARTIAL.STRUCTURE
