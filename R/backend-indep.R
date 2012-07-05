
# second prinple of CI algorithms: infer arc orientation from graph structure.
second.principle = function(x, cluster = NULL, mb, whitelist, blacklist,
    test, alpha, B = NULL, data, strict, debug = FALSE) {

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
           alpha = alpha, B = B, test = test, debug = debug))
  rownames(vs) = NULL

  if (!is.null(vs)) {

    # 3.2 sort them in p-value order.
    vs = vs[order(vs[,"max_a"], decreasing = FALSE),]
    # 3.3 apply.
    arcs = vstruct.apply(arcs = arcs, vs = vs, nodes = nodes,
             strict = strict, debug = debug)

  }#THEN

  # 4. [Remove Cycles] and 5. [Reverse Edges]
  # thou shalt not create cycles in the graph, it's acyclic!
  arcs = orient.edges(arcs = arcs, nodes = nodes, whitelist = whitelist,
           blacklist = blacklist, cluster = cluster, debug = debug,
           pass = 1)

  # 6. [Propagate Directions]
  arcs = orient.edges(arcs = arcs, nodes = nodes, whitelist = whitelist,
           blacklist = blacklist, cluster = cluster, debug = debug,
           pass = 2)

  # do a last sanity check, just in case.
  if (!is.acyclic(arcs = arcs, nodes = nodes)) {

    if (strict)
      stop("there are still cycles in the graph, aborting.")
    else
      warning("there are still cycles in the graph, fixing.")

    arcs = orient.edges(arcs = arcs, nodes = nodes, whitelist = whitelist,
             blacklist = blacklist, cluster = cluster, debug = debug,
             pass = 3)

    if (!is.acyclic(arcs = arcs, nodes = nodes))
      stop("the graph still contains cycles.")

  }#THEN

  # save the status of the learning algorithm.
  learning = list(whitelist = whitelist, blacklist = blacklist,
    test = test, args = list(alpha = alpha),
    ntests = get(".test.counter", envir = .GlobalEnv))

  # include also the number of permutations/bootstrap samples
  # if it makes sense.
  if (!is.null(B))
    learning$args$B = B

  list(learning = learning, nodes = cache.structure(nodes, arcs = arcs),
    arcs = arcs)

}#SECOND.PRINCIPLE

# build the neighbourhood of a node from the markov blanket.
neighbour = function(x, mb, data, alpha, B = NULL, whitelist, blacklist,
  backtracking = NULL, test, empty.dsep = TRUE, markov = TRUE, debug = FALSE) {

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
          function(y) { is.blacklisted(blacklist, c(x, y), both = TRUE) })]
  whitelisted = nbrhood[sapply(nbrhood,
          function(y) { is.whitelisted(whitelist, c(x, y), either = TRUE) })]

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
    cat("* learning neighbourhood of", x, ".\n")
    cat("  * blacklisted nodes: '", blacklisted, "'\n")
    cat("  * whitelisted nodes: '", whitelisted, "'\n")
    cat("  * starting with neighbourhood: '", nbrhood, "'\n")

    if (!is.null(backtracking)) {

      cat("  * known good (backtracking): '", known.good, "'.\n")
      cat("  * known bad (backtracking): '", known.bad, "'.\n")

    }#THEN

  }#THEN

  # nothing much to do, just return.
  if (length(nbrhood) <= 1)
    return(list(mb = mb[[x]], nbr = nbrhood))

  # define the backward selection heuristic.
  nbr = function(y, x, mb, test) {

    k = ifelse(empty.dsep, 0, 1)
    a = 0

    if (debug)
      cat("  * checking node", y, "for neighbourhood.\n")

    # choose the smaller set of possible d-separating nodes. 
    if (markov)
      dsep.set = smaller(mb[[x]][mb[[x]] != y], mb[[y]][mb[[y]] != x])
    else
      dsep.set = mb[[x]][mb[[x]] != y]

    if (debug)
      cat("    > dsep.set = '", dsep.set, "'\n")

    repeat {

      # create all possible subsets of the markov blanket (excluding
      # the node to be tested for exclusion) of size k.
      dsep.subsets = subsets(length(dsep.set), k, dsep.set)

      for (s in 1:nrow(dsep.subsets)) {

        if (debug)
          cat("    > trying conditioning subset '", dsep.subsets[s,], "'.\n")

        a = conditional.test(x, y, dsep.subsets[s,], data = data,
              test = test, B = B, alpha = alpha)
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

        k = k + 1

      }#THEN
      else {

        break

      }#ELSE

    }#REPEAT

  }#NBR

  # do not even try to remove whitelisted and backtracked (good) nodes.
  sapply(nbrhood[!(nbrhood %in% unique(c(whitelisted, known.good)))],
           nbr, x = x, mb = mb, test = test)

  return(list(mb = mb[[x]], nbr = nbrhood))

}#NEIGHBOUR

# detect v-structures in the graph.
vstruct.detect = function(nodes, arcs, mb, data, alpha, B = NULL, test,
    debug = FALSE) {

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

      y = tos.combs[j, 1]
      z = tos.combs[j, 2]

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

            a = conditional.test(y, z, c(dsep.subsets[s,], x), data = data,
                  test = test, B = B, alpha = alpha)
            if (debug)
              cat("    > testing", y, "vs", z, "given", c(dsep.subsets[s,], x), "(", a, ")\n")
            max_a = max(a, max_a)
            if (a > alpha) {

              if (debug)
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

    return(vs)

  }#VSTRUCT.CENTERED.ON

  sapply(nodes, vstruct.centered.on, mb = mb, data = data, simplify = FALSE)

}#VSTRUCT.DETECT

# apply v-structures to a graph.
vstruct.apply = function(arcs, vs, nodes, strict, debug = FALSE) {

  if (debug)
    cat("----------------------------------------------------------------\n")

  apply(vs, 1, function(v) {

    if (!(is.listed(arcs, v[c("y", "x")]) && is.listed(arcs, v[c("z", "x")]))) {

      if (debug) {

        cat("* not applying v-structure", v["y"], "->", v["x"], "<-", v["z"],
              "(", v["max_a"], ")\n")

      }#THEN

      if (strict)
        stop(paste("vstructure", v["y"], "->", v["x"], "<-", v["z"],
          "is not applicable, because one or both arcs are oriented",
          "in the opposite direction."))
      else
        warning(paste("vstructure", v["y"], "->", v["x"], "<-", v["z"],
          "is not applicable, because one or both arcs are oriented",
          "in the opposite direction."))

      return(NULL)

    }#THEN

    # save the updates arc set.
    temp = set.arc.direction(v["y"], v["x"], arcs)
    temp = set.arc.direction(v["z"], v["x"], temp)

    if (debug)
      cat("* applying v-structure", v["y"], "->", v["x"], "<-", v["z"],
            "(", v["max_a"], ")\n")

    assign("arcs", temp, envir = sys.frame(-2))

  })

  return(arcs)

}#VSTRUCT.APPLY

# remove arcs from the graph to make it acyclic.
orient.edges = function(arcs, nodes, whitelist, blacklist, pass, cluster,
    debug = FALSE) {

  to.be.reversed = character(0)
  narcs = nrow(arcs)
  n = 0

  # remove undirected arcs from the whitelist; their direction is
  # left to the algorithm to choose.
  if (!is.null(whitelist) && (nrow(whitelist) > 0))
    whitelist = whitelist[which.directed(whitelist),, drop = FALSE]

  # the first pass considers only directed arcs; drop the undirected ones.
  if (pass == 1) {

    which.ones = which.undirected(arcs, nodes)
    u = arcs[which.ones,, drop = FALSE]
    arcs = arcs[!which.ones,, drop = FALSE]

  }#THEN

  repeat {

    # the second pass has the main purpose to remove the propagate the directions
    # of the directed arcs; do not try to drop/reverse them again.
    if (pass == 2) {

      d = arcs[which.directed(arcs),, drop = FALSE]

      if (!is.null(whitelist))
        whitelist = arcs.rbind(whitelist, d)
      else
        whitelist = d

    }#THEN

    # cannot remove or reverse more arcs than you have.
    if (n > narcs)
      stop("too many iterations, probably would have gone on forever.")

    # check arcs between nodes which are both part of at least one cycle.
    in.cycles = is.acyclic(arcs = arcs, nodes = nodes, return.nodes = TRUE)
    to.be.checked = (arcs[, "from"] %in% in.cycles) & (arcs[, "to"] %in% in.cycles)

    if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("*", switch(pass, "detecting cycles", "propagating directions", "fixing remaining cycles"), "...\n")
      cat("  > ignored nodes:", length(nodes) - length(in.cycles), "/", length(nodes), "\n")
      cat("  > ignored arcs:", length(which(!to.be.checked)), "/", nrow(arcs), "\n")
      print(arcs[!to.be.checked, , drop = FALSE])
      cat("  > checked arcs:", length(which(to.be.checked)), "/", nrow(arcs), "\n")

    }#THEN

    # compute the cycle counter of each (relevant) arc.
    cycles = cycle.counter(arcs = arcs[to.be.checked,, drop = FALSE],
              nodes = in.cycles, cluster = cluster)

    # do not check whitelisted directed arcs.
    cycles = cycles[!which.listed(cycles[, 1:2], whitelist),, drop = FALSE]

    if (debug && (nrow(cycles) > 0)) print(cycles)

    # if there are no more arcs to check, break.
    if (nrow(cycles) == 0) break

    # if an arc belongs to a completely directed cycle but do it reverse does
    # not, drop it from the cycles data frame (this makes the second pass compliant
    # with Margaritis' algorithm specification).
    if (pass == 2) {

        amatd = arcs2amat(arcs[which.directed(arcs, nodes),, drop = FALSE], nodes)
        candidates = which.undirected(as.matrix(cycles[, 1:2, drop = FALSE]), nodes)

        for (i in which(candidates)) {

          res = !has.path(cycles[i, 2], cycles[i, 1], nodes, amatd, exclude.direct = TRUE) &&
                 has.path(cycles[i, 1], cycles[i, 2], nodes, amatd, exclude.direct = TRUE)
          candidates[i] = res

        }#FOR

        cycles = cycles[!candidates,, drop = FALSE]

    }#THEN

    if (pass != 1) {

      changed = FALSE

      # check which arcs belong to the maximum number of cycles in both
      # directions and ignore them; removing any of them does not make
      # sense.
      max.cycles = cycles[1, 3]
      which.max.cycles = (cycles[, 3] == max.cycles)
      candidates = which.undirected(as.matrix(cycles[which.max.cycles, 1:2, drop = FALSE]), nodes)

      if ((length(which(candidates)) > 0) & (nrow(cycles[!candidates,, drop = FALSE]) > 0)) {

        cycles = cycles[!candidates,, drop = FALSE]
        changed = TRUE

      }#THEN

      # if there is more than one candidate for removal, remove the one whose
      # reverse belongs to the least number of cycles.
      max.cycles = cycles[1, 3]
      which.max.cycles = (cycles[, 3] == max.cycles) & which.undirected(as.matrix(cycles[, 1:2]), nodes)

      if (length(which(which.max.cycles)) > 1) {

        # get the reverse of the arcs belonging to maximum number of cycles.
        candidates = as.matrix(cycles[which.max.cycles, 2:1, drop = FALSE])
        # get their indexes in the arc set.
        candidates = which.listed(as.matrix(cycles[, 1:2, drop = FALSE]), candidates)
        # order the arcs according to the number of cycles their reverses
        # belong to.
        cycles2 = cycles[candidates,, drop = FALSE]
        cycle.order = order(cycles2[, 3, drop = FALSE], decreasing = FALSE)
        cycles[which.max.cycles,] = cycles2[cycle.order,, drop = FALSE]
        changed = TRUE

      }#THEN

      if (debug && changed) {

        cat("  > reordered arcs (belonging to", max.cycles, "cycles each).\n")
        which.max.cycles = (cycles[, 3] == max.cycles)
        print(cycles[which.max.cycles,, drop = FALSE])

      }#THEN

    }#THEN

    # if the arc was already directed and its reverse is not blacklisted,
    # store it for future reversal.
    if (is.directed(cycles[1, c("from", "to")], arcs) &&
        !is.blacklisted(blacklist, cycles[1, c("to", "from")])) {

      if (debug)
        cat("  > arc", cycles[1, "from"], "->", cycles[1, "to"],
          "will be reversed (", cycles[1, "cycles"], " cycles ).\n")

      to.be.reversed = c(to.be.reversed, cycles[1, "from"], cycles[1, "to"])

    }#THEN
    else {

      if (debug)
        cat("  > arc", cycles[1, "from"], "->", cycles[1, "to"],
          "will be removed (", cycles[1, "cycles"], " cycles ).\n")

    }#THEN

    # remove the arc belonging to the highest number of cycles.
    arcs = arcs[!is.row.equal(arcs, cycles[1, c("from", "to")]),, drop = FALSE]

    if (debug) {

      cat("  > arcs scheduled for reversal are:\n")
      print(matrix(to.be.reversed, ncol = 2, byrow = TRUE,
              dimnames = list(NULL, c("from", "to"))))

    }#THEN

    n = n + 1

  }#REPEAT

  # add removed arcs, reversed and return.
  to.be.reversed = matrix(to.be.reversed, ncol = 2, byrow = TRUE)

  if (pass == 1)
    return(arcs.rbind(arcs.rbind(u, arcs), to.be.reversed, reverse2 = TRUE))
  else
    return(arcs.rbind(arcs, to.be.reversed, reverse2 = TRUE))

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
  counter = counter[order(counter[,"cycles"], decreasing = TRUE),, drop = FALSE]

  return(counter)

}#CYCLE.COUNTER

# emergency measures for markov blanket and neighbourhood recovery.
bn.recovery = function(bn, nodes, strict, mb = FALSE, debug = FALSE) {

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
