
# unified tabu search implementation (both optimized and by spec).
tabu.search = function(x, start, whitelist, blacklist, score,
    extra.args, max.iter, optimized, tabu, debug = FALSE) {

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
  # allocate the tabu list.
  tabu.list = vector("list", tabu)
  # maximum number of iteration the algorithm can go on without
  # improving the best network score.
  max.loss.iter = tabu
  # set the counter for suc iterations.
  loss.iter = 0
  #
  best.score = -Inf

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
    start$learning$algo = "tabu"
    start$learning$ntests = 0
    start$learning$test = score
    start$learning$args = extra.args
    start$learning$optimized = optimized

  }#THEN

  repeat {

    current = as.integer((iter - 1) %% tabu)

    # keep the best network seen so far and its score value for the
    # evaluation of the stopping rule.
    if (sum(reference.score) > best.score) {

      best.network = start
      best.score = sum(reference.score)

    }#THEN

    if (debug)
      cat("* iteration", iter, "using element", current, "of the tabu list.\n")

    # build the adjacency matrix of the current network structure.
    amat = arcs2amat(start$arcs, nodes)

    # add the hash of the network into the tabu list for future reference.
    # (BEWARE: in place modification of tabu.list!)
    .Call("tabu_hash",
          amat = amat,
          nodes = nodes,
          tabu.list = tabu.list,
          current = current,
          PACKAGE = "bnlearn")

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
    bestop = .Call("tabu_step",
                   amat = amat,
                   nodes = nodes,
                   added = to.be.added,
                   cache = cache,
                   reference = reference.score,
                   wlmat = wlmat,
                   blmat = blmat,
                   tabu.list = tabu.list,
                   current = current,
                   baseline = 0,
                   debug = debug,
                   PACKAGE = "bnlearn")

    # the value FALSE is the canary value in bestop$op meaning "no operation
    # improved the network score"; reconsider prevously discarded solutions
    # and find the one that causes the minimum decrease in the network score.
    if (bestop$op == FALSE) {

      if (loss.iter >= max.loss.iter) {

        # reset the return value to the best network ever found.
        start = best.network

        if (debug) {

          cat("----------------------------------------------------------------\n")
          cat("* maximum number of iterations without improvements reached, stopping.\n")
          cat("* best network ever seen is:\n")
          print(best.network)

        }#THEN

        break

      }#THEN
      else {

        # increase the counter of the iterations without improvements.
        loss.iter = loss.iter + 1

      }#ELSE

      if (debug) {

        cat("----------------------------------------------------------------\n")
        cat("* network score did not increase for", loss.iter,
              "iterations, looking for a minimal decrease :\n")

      }#THEN

      bestop = .Call("tabu_step",
                     amat = amat,
                     nodes = nodes,
                     added = to.be.added,
                     cache = cache,
                     reference = reference.score,
                     wlmat = wlmat,
                     blmat = blmat,
                     tabu.list = tabu.list,
                     current = current,
                     baseline = -Inf,
                     debug = debug,
                     PACKAGE = "bnlearn")

      # it might be that there are no more legal operations.
      if(bestop$op == FALSE) {

       if (debug) {
         cat("----------------------------------------------------------------\n")
         cat("* no more possible operations.\n")
         cat("@ stopping at iteration", iter, ".\n")
         }#THEN

       # reset the return value to the best network ever found.
       if (loss.iter > 0)
         start = best.network

       break

     }#THEN

    }#THEN
    else {

      if (sum(reference.score) > best.score)
        loss.iter = 0

    }#ELSE

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
      cat(sprintf("* best score up to now: %s (delta: %s)\n", format(best.score),
        format(sum(reference.score) - best.score)))

    }#THEN

    # check the current iteration index against the max.iter parameter.
    if (iter >= max.iter) {

      if (debug)
        cat("@ stopping at iteration", max.iter, ".\n")

      # reset the return value to the best network ever found.
      if (loss.iter > 0)
        start = best.network

      break

    }#THEN
    else iter = iter + 1

  }#REPEAT

  return(start)

}#TABU.SEARCH

