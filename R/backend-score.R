
# compare two network scores in an efficient way.
score.delta = function(arc, network, data, score, score.delta,
    reference.score, op, extra, debug = FALSE) {

  if (op == "reverse") {

    # do a minimal update to the network structure.
    network$nodes[[arc[1]]]$parents = c(network$nodes[[arc[1]]]$parents, arc[2])
    network$nodes[[arc[2]]]$parents =
      network$nodes[[arc[2]]]$parents[network$nodes[[arc[2]]]$parents != arc[1]]

    # compute the updated score contributions of the nodes involved.
    new.score = per.node.score(network = network, score = score,
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

    # do a minimal update to the network structure.
    if (op == "set")
      network$nodes[[arc[2]]]$parents = c(network$nodes[[arc[2]]]$parents, arc[1])
    else if (op == "drop")
      network$nodes[[arc[2]]]$parents =
        network$nodes[[arc[2]]]$parents[network$nodes[[arc[2]]]$parents != arc[1]]

    # compute the updated score contributions of arc[2].
    new.score = per.node.score(network = network, score = score,
                        nodes = arc[2], extra.args = extra, data = data)

    # update the test counter.
    assign(".test.counter", get(".test.counter", envir = .GlobalEnv) + 1,
      envir = .GlobalEnv)

    # compare the network scores.
    retval = new.score - reference.score[arc[2]]

  }#ELSE

  if (debug) {

    cat("    > delta between scores for nodes", arc, "is", retval, ".\n")

  }#THEN

  return(list(bool = (retval > score.delta), delta = retval, updates = new.score))

}#SCORE.DELTA

# create a dataframe containing the arcs to be added:
#   1               all the possibile arcs.
# - amat            exclude arcs already in the graph.
# - t(amat)         exclude the reverse of those arcs.
# - diag(n.nodes)   exclude self-loops.
# - blmat		exclude blacklisted arcs.
arcs.to.be.added = function(amat, nodes, n.nodes, blmat) {

  if (!is.null(blmat))
    amat2arcs(1 - amat - t(amat) - diag(n.nodes) - blmat, nodes)
  else
    amat2arcs(1 - amat - t(amat) - diag(n.nodes), nodes)

}#ARCS.TO.BE.ADDED

# create a dataframe containing the arcs to be dropped:
# arcs                   arcs already in the graph.
# !is.listed(whitelist)  exclude whitelisted arcs.
arcs.to.be.dropped = function(arcs, whitelist) {

  if (!is.null(whitelist))
    arcs[!which.listed(arcs, whitelist), , drop = FALSE]
  else
    arcs

}#ARCS.TO.BE.DROPPED

# create a dataframe containing the arcs to be reversed:
arcs.to.be.reversed = function(arcs, blacklist) {

  if (!is.null(blacklist))
    arcs[!which.listed(arcs[, c(2,1), drop = FALSE], blacklist), , drop = FALSE]
  else
    arcs

}#ARCS.TO.BE.REVERSED

# generic step of the standard hill-climbing algorithm.
hc.step = function(arc, op, nodes, amat, start, data, score, reference.score,
            score.delta, extra.args, debug) {

  # if the arc creates cycles, do not add/reverse it.
  if (hc.step.cycles(arc = arc, amat = amat, nodes = nodes,
        op = op, debug = debug))
          return(NA)

  if (debug) {

    if (op == "set")
      cat("  > trying to add", arc[1], "->", arc[2], ".\n")
    else if (op == "drop")
      cat("  > trying to remove", arc[1], "->", arc[2], ".\n")
    else if (op == "reverse")
      cat("  > trying to reverse", arc[1], "->", arc[2], ".\n")

  }#THEN

  # compare the scores of the two networks.
  better = score.delta(arc = arc, network = start, data = data,
             score = score, score.delta = score.delta,
             reference.score = reference.score, op = op,
             extra = extra.args, debug = debug)

  if (better$bool) {

    # store the operator and the arc for later use.
    start$lastop = c(arc[1], arc[2], op = op)

    if (debug) {

      if (op == "set") {

        cat("    @ adding", arc[1], "->", arc[2], ".\n")
        start$lastmsg = paste("adding", arc[1], "->", arc[2])

      }#THEN
      else if (op == "drop") {

        cat("    @ removing", arc[1], "->", arc[2], ".\n")
        start$lastmsg = paste("removing", arc[1], "->", arc[2])

      }#THEN
      else if (op == "reverse") {

        cat("    @ reversing", arc[1], "->", arc[2], ".\n")
        start$lastmsg = paste("reversing", arc[1], "->", arc[2])

      }#THEN

    }#THEN

    # update the network structure.
    if (op == "set")
      start$arcs = rbind(start$arcs, arc, deparse.level = 0)
    else if (op == "drop")
      start$arcs = drop.arc.backend(start$arcs, arc)
    else if (op == "reverse")
      start$arcs = reverse.arc.backend(arc[1], arc[2], start$arcs)

    start$score.delta = better$delta
    start$updates = better$updates

    return(start)

  }#THEN

  return(NA)

}#HC.STEP

# operation-aware cycle detection for the hill-climbing algorithm.
hc.step.cycles = function(arc, amat, nodes, op, debug) {

  # if the arc creates cycles, do not add/reverse it.
  if (op == "set")
    cyclic = has.path(arc[2], arc[1], nodes, amat)
  else if (op == "reverse")
    cyclic = has.path(arc[1], arc[2], nodes, amat, exclude.direct = TRUE)
  else cyclic = FALSE

  if (debug && cyclic) {

    if (op == "set")
      cat("  > trying to add", arc[1], "->", arc[2], "(cycles!).\n")
    else if (op == "reverse")
      cat("  > trying to reverse", arc[1], "->", arc[2], "(cycles!).\n")

  }#THEN

  return(cyclic)

}#HC.STEP.CYCLES

# generic step of the optimized hill-climbing algorithm.
hc.opt.step = function(arc, op, nodes, amat, start, data, score, reference.score,
    score.delta, cache, score.equivalence = FALSE, score.equivalent.cache = NULL,
    extra.args, debug) {

  # initialize the buffered score delta to NA.
  w = NA
  # initialize a flag to recall if we did hit the cache.
  cache.hit = TRUE

  # the cache is empty if the network has not yet been updated (and no
  # preseed network is used); look for a score equivalent arc whose score
  # delta is already cached.
  if (is.null(start$updates) && (op == "set")) {

    if (score.equivalence) {

      w = .Call("cache_lookup",
                arc = arc,
                cache = score.equivalent.cache,
                PACKAGE = "bnlearn")

    }#THEN

  }#THEN

  # the cache is empty if the network has not yet been updated (and no
  # preseed network is used); no need to do a lookup in this case.
  if (!is.null(start$updates)) {

    # the node(s) with the incoming arrow(s) must not have been updated
    # in the last iteration to have a valid cached score delta.
    if ((op == "reverse") && !any(arc %in% names(start$updates)) ||
        (op != "reverse") && !(arc[2] %in% names(start$updates))) {

      w =  .Call("cache_lookup",
                  arc = arc,
                  cache = cache,
                  PACKAGE = "bnlearn")

    }#THEN

    # if we are reversing an arc with with a high cached score
    # delta, we need to call score.delta() anyway beacuse we need
    # better$updates.
    if (!is.na(w) && (w > score.delta) && (op == "reverse")) w = NA

  }#THEN

  # if no cached value has been found, compute the score delta.
  if (is.na(w)) {

    if (hc.step.cycles(arc = arc, amat = amat, nodes = nodes,
          op = op, debug = debug))
            return(list(score.delta = NA))


    # reset the cache.hit flag, no cached value found.
    cache.hit = FALSE

    # compare the scores of the two networks.
    better = score.delta(arc = arc, network = start, data = data,
               score = score, score.delta = score.delta,
               reference.score = reference.score, op = op,
               extra = extra.args)

    w = better$delta

  }#THEN

  # if the arc creates cycles, do not add/reverse it.
  if (w > 0)
    if (hc.step.cycles(arc = arc, amat = amat, nodes = nodes,
          op = op, debug = debug))
            return(list(score.delta = w))

  if (debug) {

    if (op == "set")
      cat("  > trying to add", arc[1], "->", arc[2])
    else if (op == "drop")
      cat("  > trying to remove", arc[1], "->", arc[2])
    else if (op == "reverse")
      cat("  > trying to reverse", arc[1], "->", arc[2])

    cat(ifelse(cache.hit, " (cached!).\n", ".\n"))
    cat("    > delta between scores for nodes", arc, "is", w, ".\n")

  }

  # if the score delta is positive add/drop/remove the arc.
  if (w > score.delta) {

    # store the operator and the arc for later use.
    start$lastop = c(arc[1], arc[2], op = op)

    if (op == "set") {

      if (debug) {

        cat("    @ adding", arc[1], "->", arc[2], ".\n")
        start$lastmsg = paste("adding", arc[1], "->", arc[2])

      }#THEN
      start$arcs = rbind(start$arcs, arc, deparse.level = 0)
      start$updates = reference.score[arc[2]] + w
      names(start$updates) = arc[2]

    }#THEN
    else if (op == "drop") {

      if (debug)  {

        cat("    @ removing", arc[1], "->", arc[2], ".\n")
        start$lastmsg = paste("removing", arc[1], "->", arc[2])

      }#THEN
      start$arcs = drop.arc.backend(start$arcs, arc)
      start$updates = reference.score[arc[2]] + w
      names(start$updates) = arc[2]

    }#THEN
    else if (op == "reverse") {

      if (debug) {

        cat("    @ reversing", arc[1], "->", arc[2], ".\n")
        start$lastmsg = paste("reversing", arc[1], "->", arc[2])

      }#THEN
      start$arcs = reverse.arc.backend(arc[1], arc[2], start$arcs)
      start$updates = better$updates

    }#THEN

    start$score.delta = w

    return(start)

  }#ELSE
  else {

    # this is needed to update the score cache.
    return(list(score.delta = w))

  }#ELSE

}#HC.OPT.STEP

# random restart: perturb the graph and update the counter.
random.restart = function(start, x, perturb, restart, nodes, amat, score,
    extra.args, whitelist, blacklist, rebuild = FALSE, debug) {

  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* doing a random restart,", start$restart - 1, "of", restart, "left.\n")

  }#THEN

  # perturb the network.
  end = perturb.backend(network = start, iter = perturb, nodes = nodes,
          amat = amat, whitelist = whitelist, blacklist = blacklist,
          debug = debug)

  # decrement the number of random restart we still have to do.
  end$restart = end$restart - 1;

  # update the cached values of the end network.
  end$nodes = cache.structure(nodes, end$arcs)

  # update the scores of the nodes as needed.
  if (rebuild) {

    # all scores must be updated; this happens when both the network resulted
    # from the random restart is discarded and the old network is perturbed.
    end$updates = per.node.score(network = end, score = score,
                    nodes = names(end$nodes), extra.args = extra.args,
                    data = x)

  }#THEN
  else {

    # the scores whose nodes' parents changed must be updated.
    end$updates = per.node.score(network = end, score = score,
                    nodes = names(end$updates), extra.args = extra.args,
                    data = x)

  }#THEN

  return(end)

}#RANDOM.RESTART

# random restart: choose best the graph between the one found by the algorithm
# and the one found by the last random restart.
random.restart.network = function(start, restart, reference.score, debug) {

  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* current network is :\n")
    print(start)
    cat("* current score:", sum(reference.score), "\n")

  }#THEN

  if (start$restart == restart) {

    # this condition is triggered on the first random restart.
    # save the (hopefully) optimal model found by the algorithm, which
    # is used as a starting network for the random restarts, along with
    # its score.
    start$learning$nodes = start$nodes
    start$learning$arcs = start$arcs
    start$learning$score = sum(reference.score)

  }#THEN
  else {

    # if the network found by the algorithm is better, use that one as
    # starting network for the next random restart; use the old one
    # otherwise.
    if (start$learning$score > sum(reference.score)) {

      start$arcs = start$learning$arcs
      start$nodes = start$learning$nodes

      if (debug) {

        cat("* the old network was better, discarding the current one.\n")
        cat("* now the current network is :\n")
        print(start)
        cat("* current score:", start$learning$score, "\n")

       }#THEN

    }#THEN
    else {

      start$learning$nodes = start$nodes
      start$learning$arcs = start$arcs
      start$learning$score = sum(reference.score)

    }#ELSE

  }#ELSE

  return(start)

}#RANDOM.RESTART.NETWORK

