
# second prinple of CI algorithms: infer arc orientation from graph structure.
second.principle = function(x, cluster = NULL, mb, nodes, whitelist, blacklist, 
  test, alpha, data, strict, direction, debug) {

  # hope it's never called ...
  mb = nbr.recovery(mb, nodes = nodes, strict = strict, debug = debug)

  # build a list of the undirected arcs in the graph, using the neighbourhoods
  # detected in markov.blanket().
  arcs = mb2arcs(mb, nodes) 

  # apply blacklist to the arc list.
  arcs = arcs[!apply(arcs, 1, function(x){ is.blacklisted(blacklist, x) }),]

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
  # thou shalt not create loops in the graph, it's acyclic!
  arcs = orient.edges(arcs = arcs, nodes = nodes, 
           whitelist = whitelist, debug = debug)

  # arcs whitelisted in one direction (i.e. a -> b but not b -> a) are
  # never checked for loops, because they are definitely not
  # undirected.
  # do the check now.
  if (!is.acyclic(list(nodes = mb, arcs = arcs)))
    stop("the graph contains cycles because of whitelisted nodes.")

  # 6. [Propagate Directions]
  arcs = propagate.directions(arcs = arcs, nodes = nodes, debug = debug)

  # EXTRA [ESP]
  if (direction)
    arcs = set.directions(arcs = arcs, data = x, test = test, 
             alpha = alpha, cluster = cluster, debug = debug)

  list(nodes = mb, arcs = arcs, whitelist = whitelist, 
    blacklist = blacklist, test = test, alpha = alpha, 
    ntests = get("test.counter", envir = .GlobalEnv))

}#SECOND.PRINCIPLE

propagate.directions = function(arcs, nodes, debug) {

  undirected.arcs = arcs[which.undirected(arcs),]

  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* propagating directions for the following undirected arcs:\n")

  }#THEN

  if (nrow(undirected.arcs) == 0) {

    if (debug) cat("  > nothing to do.\n")
 
  }#THEN
  else apply(undirected.arcs, 1, function(arc) {

    if (has.path(arcs, nodes, arc[1], arc[2])) {

      assign("arcs", arcs[!is.row.equal(arcs, arc[c(2,1)]),],
        envir = sys.frame(-2))

      if (debug)
        cat("  > there's a path from", arc[1], "to", arc[2], ".\n")
        cat("    > removing", arc[2], "->", arc[1], ".\n")

    }#THEN

    if (debug) 
      cat("  > no path from", arc[1], "to", arc[2], ".\n")

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

          a = conditional.test(x, y, dsep.subsets[s,], data = data, test = test)
          if (a > alpha) {

            if (debug) {

              cat("    > node", y, "is not a neighbour of", x, ". ( p-value:", a, ")\n")
              cat("    > conditioning subset: '", dsep.subsets[s,], "'\n")

            }#THEN

            # update the neighbourhood.
            assign("nbrhood", nbrhood[nbrhood != y], envir = sys.frame(-3) )

            break

          }#THEN
          else if (debug) {

            cat("    > trying conditioning subset '", dsep.subsets[s,], 
                "' ( p-value:", a, ")\n")

          }#ELSE

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

      assign("arcs", arcs[!is.row.equal(arcs, v[c("x", "y")]), ],
        envir = sys.frame(-2))
      assign("arcs", arcs[!is.row.equal(arcs, v[c("x", "z")]), ],
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
orient.edges = function(arcs, nodes, whitelist, debug) {

  to.be.reversed = matrix(rep(0, 2), ncol = 2, 
                     dimnames = list(c(), c("from", "to")))[0,]
  narcs = nrow(arcs)
  n = 0

  # remove undirected arcs from the whitelist; their direction is
  # left to the algorithm to choose.
  if (!is.null(whitelist) && (nrow(whitelist) > 0))
    whitelist = whitelist[!which.undirected(whitelist), ]

  # be really sure that the whitelist is still a amtrix.
  if (is.character(whitelist))
    whitelist = matrix(whitelist, ncol = 2, dimnames = list(c(), c("from", "to")))

  repeat {

    # cannot remove or reverse more arcs than you have.
    if (n > narcs)
      stop("too many iterations, probably would have gone on forever.")

    if (debug) {
      
      cat("----------------------------------------------------------------\n")
      cat("* detecting loops ...\n")

    }#THEN

    loops = loop.counter(arcs, nodes)

    # do not check arcs whitelisted directed arcs.
    loops = loops[!which.listed(loops[, 1:2], whitelist), ]

    if (debug) print(loops)

    if (all(loops[,3] == 0)) break

    arcs = arcs[!is.row.equal(arcs, loops[1, c("from", "to")]),]

    # if the arc was already directed, store it for future reversal.
    if (is.listed(arcs, loops[1, c("from", "to")], both = TRUE))
      to.be.reversed = rbind(to.be.reversed, loops[1, c("from", "to")])          

    if (debug) {

      cat("  > arc", loops[1, "from"], "->", loops[1, "to"], 
        "will be removed (", loops[1, "loops"], " cycles ).\n")

      cat("  > arcs scheduled for reversal are:\n")
      print(to.be.reversed)

    }#THEN

    n = n + 1

  }#REPEAT

  # add removed arcs, reversed.
  arcs = rbind(arcs, to.be.reversed[, c(2,1)])

  arcs

}#ORIENT.EDGES

# count how many loops each arc is part of.
loop.counter = function(arcs, nodes) {

  # build the adjacency matrix.
  amat = arcs2amat(arcs, nodes)
  # loop counter
  counter = cbind(arcs, loops = rep(0, nrow(arcs)))

  apply(arcs, 1,

    function(arc) {

      buffer = arc
      dim(buffer) = c(1,2)

      for (i in 1:(length(nodes) - 1)) {

        buffer = apply(buffer, 1, function(buf) {

          next.buffer = NULL

          # get the next nodes in this path. NA means no more nodes, i.e no
          # outgoing arc from the last node of the path.
          if (!is.na(buf[length(buf)]))
            next.one = names(which(amat[buf[length(buf)],] == 1))
          else
            next.one = NA

          # discard backward steps due to unoriented arcs, i.e. a -> b -> a
          next.one = next.one[next.one != buf[length(buf) - 1]]

          # if there's noone left, set next.one to NA.
          if (length(next.one) == 0) next.one = NA

          # rebuild the buffer.
          if (is.null(next.buffer)) 
            next.buffer = t(sapply(next.one, function(x) { c(buf, x) } ))
          else 
            next.buffer = rbind(sapply(next.one, function(x) { c(buf, x) } ), next.buffer)

          next.buffer

        })

        # sanity check: sometimes something goes wrong and the buffer is a list.
        if (is.list(buffer)) buffer = t(do.call("rbind", unique(buffer)))
        # hey, it's a matrix!
        buffer = matrix(t(buffer), ncol = 2 + i)

      }#FOR

      # find out how many times arc["from"] is present in any path
      loops = apply(buffer, 1, function(path) {length(which(path == arc["from"]))})

      # set the loop counter.
      counter[is.row.equal(counter[, 1:2], arc), "loops"] = length(which(loops > 1))
      assign("counter", counter, envir=sys.frame(-2))

    })

    # return the sorted loop counter.
    counter[order(as.numeric(counter[,"loops"]), decreasing = TRUE),]

}#LOOP.COUNTER

set.directions = function(arcs, data, test, alpha, cluster, debug) {

  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* infer arc direction.\n")

  }#THEN

  # get the undirected arcs.
  undirected.arcs = arcs[which.undirected(arcs), ]

  if (is.null(cluster)) {

    a = apply(undirected.arcs, 1, function(arc) {

      # you can't help but notice nodes connecetd by undirected arcs are
      # included, too? wondwer why?
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
      # included, too? wondwer why?
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

    # if both arcs have the same p-value, do not remove them from the arc list.
    if (da.best[3] != lookup.table[is.row.equal(lookup.table[, 1:2], da.best[c(2, 1)]), 3]) {

      # remove its reverse from the arc list.
      arcs = arcs[!is.row.equal(arcs, da.best[c(2:1)]), ]

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

mb.recovery = function(mb, nodes, strict, debug) { 

  if (is.symmetric(mb2amat(mb, redux = TRUE)))
    return(mb)
 
  # if the learning correctness is relaxed, do not raise an error 
  # if the markov blanket are not symmetric but the neighbourhoods 
  # are. 
  # NOTE: backtracking guarantees symmetry. 
  if (strict) { 

    if (debug) {

      m = mb2amat(mb, redux = TRUE)
      rownames(m) = colnames(m)

      cat("----------------------------------------------------------------\n") 
      cat("> markov blankets as an adjacency matrix:\n")
      print(m)
      cat("----------------------------------------------------------------\n") 
      cat("> markov blankets as a list:\n")
      print(mb) 

    }#THEN
 
    stop("markov blankets are not symmetric.") 
 
  }#THEN 
  else { 
 
    # fix the markov blankets: 
    # if X \not\in MB(Y) => Y \not\in MB(X) 
    mb = apply(mb2amat(mb, redux = TRUE) *  t(mb2amat(mb, redux = TRUE)), 1, 
            function(x) { nodes[as.logical(x)] }) 
    names(mb) = nodes 
 
    if (debug) { 
  
      cat("----------------------------------------------------------------\n") 
      cat("* fixing symmetry in the markov blankets.\n") 
      print(mb) 
 
    }#THEN 
 
    warning("markov blankets are not symmetric.") 
 
  }#ELSE 
 
  mb 
 
}#MB.RECOVERY 

nbr.recovery = function(mb, nodes, strict, debug) { 

  # you should not be here! no, really!
  if (is.symmetric(nbr2amat(mb)))
    return(mb)

  if (strict) {

    if (debug) {

      m = nbr2amat(mb)
      rownames(m) = colnames(m)

      cat("----------------------------------------------------------------\n")
      cat("> neighbourhoods as an adjacency matrix:\n")
      print(m)
      cat("----------------------------------------------------------------\n")
      cat("> neighbourhoodss as a list:\n")
      print(mb)

     }#THEN

     stop("neighbourhoods are not symmetric.")

  }#THEN
  else {

    # if X \not\in NBR(Y) => Y \not\in NBR(X) 
    nbrhood = apply(nbr2amat(mb) *  t(nbr2amat(mb)), 1,
                function(x) { nodes[as.logical(x)] })
    names(nbrhood) = nodes

    for (n in nodes)
      mb[[n]]$nbr = nbrhood[[n]]

    if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("* fixing symmetry in the neighbourhoods.\n")
      print(mb)

    }#THEN 

    warning("neighbourhoods are not symmetric.")

  }#ELSE

  mb

}#NBR.RECOVERY
