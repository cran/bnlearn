
# unified hill climbing implementation (both optimized and by spec).
hill.climbing = function(x, start, whitelist, blacklist, score,
    extra.args, restart, perturb, max.iter, optimized, debug = FALSE) {

  # cache nodes' labels.
  nodes = names(x)
  # cache the number of nodes.
  n.nodes = length(nodes)
  # set the iteration counter.
  iter = 1
  # check whether the score is score-equivalent.
  score.equivalence = score %in% score.equivalent.scores
  # allocate the cache matrix.
  cache = matrix(0, nrow = n.nodes, ncol = n.nodes)
  # nodes to be updated (all of them in the first iteration).
  updated = seq_len(n.nodes) - 1L
  # set the number of random restarts.
  restart.counter = restart

  # set the reference score.
  reference.score = per.node.score(network = start, score = score,
                      nodes = nodes, extra.args = extra.args, data = x)

  # convert the blacklist to an adjacency matrix for easy use.
  if (!is.null(blacklist))
    blmat = arcs2amat(blacklist, nodes)
  else
    blmat = matrix(0L, nrow = n.nodes, ncol = n.nodes)

  # convert the whitelist to an adjacency matrix for easy use.
  if (!is.null(whitelist))
    wlmat = arcs2amat(whitelist, nodes)
  else
    wlmat = matrix(0L, nrow = n.nodes, ncol = n.nodes)

  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* starting from the following network:\n")
    print(start)
    cat("* current score:", sum(reference.score), "\n")
    cat("* whitelisted arcs are:\n")
    if (!is.null(whitelist)) print(whitelist)
    cat("* blacklisted arcs are:\n")
    if (!is.null(blacklist)) print(blacklist)

    # set the metadata of the network; othewise the debugging output is
    # confusing and not nearly as informative.
    start$learning$algo = "hc"
    start$learning$ntests = 0
    start$learning$test = score
    start$learning$args = extra.args
    start$learning$optimized = optimized

  }#THEN

  repeat {

    # build the adjacency matrix of the current network structure.
    amat = arcs2amat(start$arcs, nodes)

    # set up the score cache (BEWARE: in place modification!).
    .Call("hc_cache_fill",
          nodes = nodes,
          data = x,
          network = start,
          score = score,
          extra = extra.args,
          reference = reference.score,
          equivalence = score.equivalence && optimized,
          updated = (if (optimized) updated else seq(length(nodes)) - 1L),
          env = environment(),
          amat = amat,
          cache = cache,
          blmat = blmat,
          debug = debug,
          PACKAGE = "bnlearn")

    # select which arcs should be tested for inclusion in the graph (hybrid
    # learning algorithms should hook the restrict phase here).
    to.be.added = arcs.to.be.added(amat = amat, nodes = nodes,
                    blacklist = blmat, whitelist = NULL, arcs = FALSE)

    # get the best arc addition/removal/reversal.
    bestop = .Call("hc_opt_step",
                   amat = amat,
                   nodes = nodes,
                   added = to.be.added,
                   cache = cache,
                   reference = reference.score,
                   wlmat = wlmat,
                   blmat = blmat,
                   debug = debug,
                   PACKAGE = "bnlearn")

    # the value FALSE is the canary value in bestop$op meaning "no operation
    # improved the network score"; break the loop.
    if (bestop$op == FALSE) {

      if ((restart.counter >= 0) && (restart > 0)) {

        if (restart.counter == restart) {

          # store away the learned network at the first restart.
          best.network = start
          best.network.score = sum(reference.score)

        }#THEN
        else {

          # if the network found by the algorithm is better, use that one as
          # starting network for the next random restart; use the old one
          # otherwise.
          if (best.network.score > sum(reference.score)) {

            start = best.network

            if (debug) {

              cat("----------------------------------------------------------------\n")
              cat("* the old network was better, discarding the current one.\n")
              cat("* now the current network is :\n")
              print(start)
              cat("* current score:", best.network.score, "\n")

            }#THEN

          }#THEN
          else {

            # store away the learned network at the first restart.
            best.network = start
            best.network.score = sum(reference.score)

          }#ELSE

        }#ELSE

        # this is the end; if the network found by the algorithm is the best one
        # it's now the 'end' one once more.
        if (restart.counter == 0) break

        # don't try to do anything if there are no more iterations left.
        if (iter >= max.iter) {

          if (debug)
            cat("@ stopping at iteration", max.iter, "ignoring random restart.\n")

          break

        }#THEN
        # increment the iteration counter.
        iter = iter + 1

        # decrement the number of random restart we still have to do.
        restart.counter = restart.counter - 1

        if (debug) {

          cat("----------------------------------------------------------------\n")
          cat("* doing a random restart,", restart.counter, "of", restart, "left.\n")

        }#THEN

        # perturb the network.
        start = perturb.backend(network = start, iter = perturb, nodes = nodes,
                amat = amat, whitelist = whitelist, blacklist = blacklist,
                debug = debug)

        # update the cached values of the end network.
        start$nodes = cache.structure(nodes, arcs = start$arcs)

        # update the scores of the nodes as needed.
        if (best.network.score > sum(reference.score)) {

          # all scores must be updated; this happens when both the network
          # resulted from the random restart is discarded and the old network
          # is perturbed.
          reference.score = per.node.score(network = start, score = score,
                          nodes = nodes, extra.args = extra.args, data = x)

          updated = which(nodes %in% nodes) - 1L

        }#THEN
        else {

          # the scores whose nodes' parents changed must be updated.
          reference.score[names(start$updates)] =
            per.node.score(network = start, score = score, nodes = names(start$updates),
               extra.args = extra.args, data = x)

          updated = which(nodes %in% names(start$updates)) - 1L

        }#THEN

        # nuke start$updates from orbit.
        start$updates = NULL

        next

      }#THEN
      else
        break

    }#THEN

    # update the network structure.
    start = arc.operations(start, from = bestop$from, to = bestop$to,
              op = bestop$op, check.cycles = FALSE, update = TRUE,
              debug = FALSE)

    # set the nodes whose cached score deltas are to be updated.
    if (bestop$op == "reverse")
      updated = which(nodes %in% c(bestop$from, bestop$to)) - 1L
    else
      updated = which(nodes %in% bestop$to) - 1L

    if (debug) {

      # update the test counter of the network; very useful to check how many
      # score comparison has been done up to now.
      start$learning$ntests = get(".test.counter", envir = .GlobalEnv)

      cat("----------------------------------------------------------------\n")
      cat("* best operation was: ")
      if (bestop$op == "set")
        cat("adding", bestop$from, "->", bestop$to, ".\n")
      else if (bestop$op == "drop")
        cat("removing", bestop$from, "->", bestop$to, ".\n")
      else
        cat("reversing", bestop$from, "->", bestop$to, ".\n")
      cat("* current network is :\n")
      print(start)
      cat("* current score:", sum(reference.score), "\n")

    }#THEN

    # check the current iteration index against the max.iter parameter.
    if (iter >= max.iter) {

      if (debug)
        cat("@ stopping at iteration", max.iter, ".\n")

      break

    }#THEN
    else iter = iter + 1

  }#REPEAT

  return(start)

}#HILL.CLIMBING

