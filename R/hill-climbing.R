
hill.climbing = function(x, start, whitelist, blacklist, score,
    extra.args, restart, perturb, max.iter, debug) {

  # cache nodes' labels.
  nodes = names(x)
  # cache the number of nodes.
  n.nodes = length(nodes)
  # set the number of random restarts.
  start$restart = restart
  # create the prospective model.
  end = start
  # set the score delta.
  end$score.delta = 0
  # set the iteration counter.
  iter = 1

  # wrapper function for the generic hill-climbing step.
  step = function(arc, op) {

    start = hc.step(arc = arc, nodes = nodes, amat = amat,
              start = start, data = x, score = score,
              reference.score = reference.score, op = op,
              score.delta = end$score.delta,
              extra.args = extra.args, debug = debug)

    if (is(start, "bn"))
      assign("end", start, envir = sys.frame(-2))

  }#STEP

  # set the reference score.
  reference.score = per.node.score(network = start, score = score,
                      nodes = nodes, extra.args = extra.args, data = x)

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

  }#THEN

  if (!is.null(blacklist))
    blmat = arcs2amat(blacklist, nodes)
  else
    blmat = NULL

  repeat {

    amat = arcs2amat(start$arcs, nodes)

    to.be.added = arcs.to.be.added(amat = amat, nodes = nodes, n.nodes = n.nodes,
                    blmat = blmat)

    if (nrow(to.be.added) > 0) {

      if (debug) {

        cat("----------------------------------------------------------------\n")
        cat("* trying to add one of", nrow(to.be.added), "arcs.\n")

      }#THEN

      apply(to.be.added, 1, step, op = "set")

    }#THEN
    else if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("* no arc to be added.\n")

    }#THEN

    if (nrow(start$arcs) > 0) {

      to.be.dropped = arcs.to.be.dropped(arcs = start$arcs, whitelist = whitelist)

      # if there is any arc in the graph, try to remove it.
      if (nrow(to.be.dropped) > 0) {

        if (debug) {

          cat("----------------------------------------------------------------\n")
          cat("* trying to remove one of", nrow(start$arcs), "arcs.\n")

        }#THEN

        apply(to.be.dropped, 1, step, op = "drop")

      }#THEN
      else if (debug) {

        cat("----------------------------------------------------------------\n")
        cat("* no arc to be removed.\n")

      }#THEN

      to.be.reversed = arcs.to.be.reversed(arcs = start$arcs, blacklist = blacklist)

      if (nrow(to.be.reversed) > 0) {

        if (debug) {

          cat("----------------------------------------------------------------\n")
          cat("* trying to reverse one of", nrow(start$arcs), "arcs.\n")

        }#THEN

        apply(to.be.reversed, 1, step, op = "reverse")

      }#THEN
      else if (debug) {

        cat("----------------------------------------------------------------\n")
        cat("* no arc to be reversed.\n")

      }#THEN

    }#THEN
    else if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("* no arc to be removed or reversed.\n")

    }#THEN

    # if the network is the same as the last iteration, give up.
    if (identical(start$arcs, end$arcs)) {

      # the likelihood score has this nasty habit to be too small for its own
      # good, and is often rounded to zero due to machine precision limit. If
      # that's the case, switch to the loglikelihood.
      if ((score == "lik") && (sum(reference.score) == 0)) {

        warning("switching to loglikelihood due to machine precision limits.")

        score = "loglik"
        reference.score = sapply(names(start$nodes), loglik.node,
                            x = start, data = x)

      }#THEN
      else if ((start$restart >= 0) && (restart > 0)) {

        start = random.restart.network(start = start, restart = restart,
                  reference.score = reference.score, debug = debug)

        # this is the end; if the network found by the algorithm is the best one
        # it's now the 'end' one once more.
        if (start$restart == 0) break

        # don't try to do anything if there are no more iterations left.
        if (iter == max.iter) {

          if (debug)
            cat("@ stopping at iteration", max.iter, "ignoring random restart.\n")

          break

        }#THEN

        # do the random restart.
        end = random.restart(start = start, x = x, restart = restart,
                perturb = perturb, nodes = nodes, amat = amat, score = score,
                extra.args = extra.args, debug = debug,
                whitelist = whitelist, blacklist = blacklist,
                rebuild = (start$learning$score > sum(reference.score)))

      }#THEN
      else
        break

    }#THEN

    # update the cached values of the end network.
    end$nodes = cache.structure(nodes, arcs = end$arcs)
    # reset the score delta.
    end$score.delta = 0
    # update the starting network for the next iteration.
    start = end
    # reset the reference score.
    if (!is.null(end$updates))
      reference.score[names(end$updates)] = end$updates

    if (debug) {

      # update the test counter of the network; very useful to check how many
      # score comparison has been done up to now.
      start$learning$ntests = get(".test.counter", envir = .GlobalEnv)

      cat("----------------------------------------------------------------\n")
      cat("* best operation was:", start$lastmsg, ".\n")
      cat("* current network is :\n")
      print(start)
      cat("* current score:", sum(reference.score), "\n")

    }#THEN

    # check the current iteration index against the max.iter parameter.
    if (iter == max.iter) {

      if (debug)
        cat("@ stopping at iteration", max.iter, ".\n")

      break

    }#THEN
    else iter = iter + 1

  }#REPEAT

  # remove all the extra elements from the return value.
  end$updates = end$score.delta = end$learning$score = end$restart = NULL

  end

}#HILL.CLIMBING

hill.climbing.optimized = function(x, start, whitelist, blacklist, score,
    extra.args, restart, perturb, max.iter, debug) {

  # cache nodes' labels.
  nodes = names(x)
  # cache the number of nodes.
  n.nodes = length(nodes)
  # set the number of random restarts.
  start$restart = restart
  # create the prospective model.
  end = start
  # set the score delta.
  end$score.delta = 0
  # set the iteration counter.
  iter = 1
  # check whether the score is score-equivalent.
  score.equivalence = score %in% score.equivalent.scores
  # set up the score caches.
  score.cache.add = score.cache.drop = score.cache.reverse =
    data.frame(from = character(0), to = character(0),
      delta = numeric(0), stringsAsFactors = FALSE)

  if (score.equivalence) {

    # set up the score-equivalence 'shadow' cache.
    shadow.cache.add = data.frame(from = character(0), to = character(0),
                         delta = numeric(0), stringsAsFactors = FALSE)

  }#THEN

  # function to add an arc to the graph.
  step = function(arc, op, cache) {

    start = hc.opt.step(arc = arc, nodes = nodes, amat = amat,
              start = start, data = x, score = score,
              reference.score = reference.score, op = op,
              score.delta = end$score.delta, cache = cache,
              score.equivalence = score.equivalence,
              score.equivalent.cache = shadow.cache.add,
              extra.args = extra.args, debug = debug)

    if (is(start, "bn"))
      assign("end", start, envir = sys.frame(-3))

    # update the score-equivalence 'shadow' cache until I have a real cache.
    if ((nrow(score.cache.add) == 0) && score.equivalence && (op == "set")) {

      assign("shadow.cache.add", rbind(shadow.cache.add, data.frame(from = arc[2],
        to = arc[1], delta = start$score.delta, stringsAsFactors = FALSE)),
        envir = sys.frame(-3))

    }#THEN

    return(start$score.delta)

  }#STEP

  # set the reference score.
  reference.score = per.node.score(network = start, score = score,
                      nodes = nodes, extra.args = extra.args, data = x)

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

  }#THEN

  if (!is.null(blacklist))
    blmat = arcs2amat(blacklist, nodes)
  else
    blmat = NULL

  repeat {

    amat = arcs2amat(start$arcs, nodes)

    to.be.added = arcs.to.be.added(amat = amat, nodes = nodes, n.nodes = n.nodes,
                    blmat = blmat)

    if (nrow(to.be.added) > 0) {

      if (debug) {

        cat("----------------------------------------------------------------\n")
        cat("* trying to add one of", nrow(to.be.added), "arcs.\n")

      }#THEN

      # try to add any available arc.
      score.cache.add = data.frame(to.be.added,
          delta = apply(to.be.added, 1, step, op = "set",
          cache = score.cache.add), stringsAsFactors = FALSE)

    }#THEN
    else if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("* no arc to be added.\n")

    }#THEN

    if (nrow(start$arcs) > 0) {

      to.be.dropped = arcs.to.be.dropped(arcs = start$arcs, whitelist = whitelist)

      # if there is any arc in the graph, try to remove it.
      if (nrow(to.be.dropped) > 0) {

        if (debug) {

          cat("----------------------------------------------------------------\n")
          cat("* trying to remove one of", nrow(start$arcs), "arcs.\n")

        }#THEN

        # try to remove any arc in the graph.
        score.cache.drop = data.frame(to.be.dropped,
            delta = apply(to.be.dropped, 1, step, op = "drop",
            cache = score.cache.drop), stringsAsFactors = FALSE)

      }#THEN
      else if (debug) {

        cat("----------------------------------------------------------------\n")
        cat("* no arc to be removed.\n")

      }#THEN

      to.be.reversed = arcs.to.be.reversed(arcs = start$arcs, blacklist = blacklist)

      if (nrow(to.be.reversed) > 0) {

        if (debug) {

          cat("----------------------------------------------------------------\n")
          cat("* trying to reverse one of", nrow(start$arcs), "arcs.\n")

        }#THEN

        # try to remove any arc in the graph.
        score.cache.reverse = data.frame(to.be.reversed,
            delta = apply(to.be.reversed, 1, step, op = "reverse",
            cache = score.cache.reverse), stringsAsFactors = FALSE)

      }#THEN
      else if (debug) {

        cat("----------------------------------------------------------------\n")
        cat("* no arc to be reversed.\n")

      }#THEN

    }#THEN
    else if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("* no arc to be removed or reversed.\n")

    }#THEN

    # if the network is the same as the last iteration, give up.
    if (identical(start$arcs, end$arcs)) {

      # the likelihood score has this nasty habit to be too small for its own
      # good, and is often rounded to zero due to machine precision limit. If
      # that's the case, switch to the loglikelihood.
      if ((score == "lik") && (sum(reference.score) == 0)) {

        warning("switching to loglikelihood due to machine precision limits.")

        score = "loglik"
        reference.score = sapply(names(start$nodes), loglik.node,
                            x = start, data = x)

      }#THEN
      else if ((start$restart >= 0) && (restart > 0)) {

        start = random.restart.network(start = start, restart = restart,
                  reference.score = reference.score, debug = debug)

        # this is the end; if the network found by the algorithm is the best one
        # it's now the 'end' one once more.
        if (start$restart == 0) break

        # don't try to do anything if there are no more iterations left.
        if (iter == max.iter) {

          if (debug)
            cat("@ stopping at iteration", max.iter, "ignoring random restart.\n")

          break

        }#THEN

        # do the random restart.
        end = random.restart(start = start, x = x, restart = restart,
                perturb = perturb, nodes = nodes, amat = amat, score = score,
                extra.args = extra.args, debug = debug,
                whitelist = whitelist, blacklist = blacklist,
                rebuild = (start$learning$score > sum(reference.score)))

      }#THEN
      else
        break

    }#THEN

    # update the cached values of the end network.
    end$nodes = cache.structure(nodes, arcs = end$arcs)
    # reset the score delta.
    end$score.delta = 0
    # update the starting network for the next iteration.
    start = end
    # reset the reference score.
    if (!is.null(end$updates))
      reference.score[names(end$updates)] = end$updates

    if (debug) {

      # update the test counter of the network; very useful to check how many
      # score comparison has been done up to now.
      start$learning$ntests = get(".test.counter", envir = .GlobalEnv)

      cat("----------------------------------------------------------------\n")
      cat("* best operation was:", start$lastmsg, ".\n")
      cat("* current network is :\n")
      print(start)
      cat("* current score:", sum(reference.score), "\n")

    }#THEN

    # check the current iteration index against the max.iter parameter.
    if (iter == max.iter) {

      if (debug)
        cat("@ stopping at iteration", max.iter, ".\n")

      break

    }#THEN
    else iter = iter + 1

  }#REPEAT

  # remove all the extra elements from the return value.
  end$updates = end$score.delta = end$learning$score = end$restart = NULL

  end

}#HILL.CLIMBING.OPTIMIZED

