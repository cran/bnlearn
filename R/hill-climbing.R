
hill.climbing = function(x, start, whitelist, blacklist, score,
    extra.args, debug) {

  # cache nodes' labels.
  nodes = colnames(x)
  # cache the number of nodes.
  n.nodes = length(nodes)
  # create the prospective model.
  end = start
  # set the score delta.
  end$score.delta = 0

  # function to add an arc to the graph.
  add.step = function(arc) {

    # if the arc creates cycles, do not add it.
    if (has.path(arc[2], arc[1], nodes, amat)) {

      if (debug) cat("  > trying to add", arc[1], "->", arc[2], "(cycles!).\n")

      return(NA)

    }#THEN

    if (debug)
      cat("  > trying to add", arc[1], "->", arc[2], ".\n")

    # compare the scores of the two networks.
    better = score.delta(arc = arc, network = start, data = x,
               score = score, score.delta = end$score.delta,
               reference.score = reference.score, op = "set",
               extra = extra.args, debug = debug)

    if (better$bool) {

      if (debug) cat("    @ adding", arc[1], "->", arc[2], ".\n")

      # no need to call set.arc.direction, the arc is not present.
      start$arcs = rbind(start$arcs, arc, deparse.level = 0)
      start$score.delta = better$delta
      start$updates = better$updates

      assign("end", start, envir = sys.frame(-2))

    }#THEN

  }#ADD.STEP

  # function to drop an arc from the graph.
  drop.step = function(arc) {

    if (debug)
      cat("  > trying to remove", arc[1], "->", arc[2], ".\n")

    # compare the scores of the two networks.
    better = score.delta(arc = arc, network = start, data = x,
               score = score, score.delta = end$score.delta,
               reference.score = reference.score, op = "drop",
               extra = extra.args, debug = debug)

    if (better$bool) {

      if (debug) cat("    @ removing", arc[1], "->", arc[2], ".\n")

      # drop the arc from the network.
      start$arcs = drop.arc.backend(start$arcs, arc)
      start$score.delta = better$delta
      start$updates = better$updates

      assign("end", start, envir = sys.frame(-2))

    }#THEN

  }#DROP.STEP

  # function to reverse an arc in the graph.
  reverse.step = function(arc) {

    # if the arc creates cycles, do not reverse it.
    if (has.path(arc[1], arc[2], nodes, amat, exclude.direct = TRUE)) {

      if (debug) cat("  > trying to reverse", arc[1], "->", arc[2], "(cycles!).\n")

      return(FALSE)

    }#THEN

    if (debug)
      cat("  > trying to reverse", arc[1], "->", arc[2], ".\n")

    # compare the scores of the two networks.
    better = score.delta(arc = arc, network = start, data = x,
               score = score, score.delta = end$score.delta,
               reference.score = reference.score, op = "reverse",
               extra = extra.args, debug = debug)

    if (better$bool) {

      if (debug) cat("    @ reversing", arc[1], "->", arc[2], ".\n")

      start$arcs = reverse.arc.backend(arc[1], arc[2], start$arcs)
      start$score.delta = better$delta
      start$updates = better$updates

      assign("end", start, envir = sys.frame(-2))

    }#THEN

    return(better$delta)

  }#REVERSE.STEP

  # set the reference score.
  if (score == "k2") {

    reference.score = sapply(names(start$nodes), dirichlet.node,
                        x = start, data = x)

  }#THEN
  else if (score %in% c("bde", "dir")) {

    reference.score = sapply(names(start$nodes), dirichlet.node,
                        x = start, imaginary.sample.size = extra.args$iss,
                        data = x)

  }#THEN
  else if (score %in% c("lik", "loglik")) {

    reference.score = sapply(names(start$nodes), loglik.node,
                        x = start, data = x)
    if (score == "lik") reference.score = exp(reference.score)

  }#THEN
  else if (score %in% c("aic", "bic")) {

    reference.score = sapply(names(start$nodes), aic.node, x = start, data = x,
                        k = extra.args$k)

  }#THEN
  else if (score == "bge") {

    reference.score = sapply(names(start$nodes), bge.node,
                        x = start, imaginary.sample.size = extra.args$iss,
                        phi = extra.args$phi, data = x)

  }#THEN

  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* starting from the following network:\n")
    print(start)
    cat("* current score:", sum(reference.score), "\n")
    cat("* whitelisted arcs are:\n")
    if (!is.null(whitelist)) print(whitelist)
    cat("* blacklisted arcs are:\n")
    if (!is.null(blacklist)) print(blacklist)

  }#THEN

  if (!is.null(blacklist))
    blmat = arcs2amat(blacklist, nodes)

  repeat {

    amat = arcs2amat(start$arcs, nodes)

    # create the list of arcs to be added.
    #   1               all the possibile arcs.
    # - amat            exclude arcs already in the graph.
    # - t(amat)         exclude the reverse of those arcs.
    # - diag(n.nodes)   exclude self-loops.
    # - blmat		exclude blacklisted arcs.
    if (!is.null(blacklist))
      to.be.added = amat2arcs(1 - amat - t(amat) - diag(n.nodes) - blmat, nodes)
    else
      to.be.added = amat2arcs(1 - amat - t(amat) - diag(n.nodes), nodes)

    if (nrow(to.be.added) > 0) {

      if (debug) {

        cat("----------------------------------------------------------------\n")
        cat("* trying to add one of", nrow(to.be.added), "arcs.\n")

      }#THEN

      apply(to.be.added, 1, add.step)

    }#THEN
    else if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("* no arc to be added.\n")

    }#THEN

    if (nrow(start$arcs) > 0) {

      # create a list of arcs to be dropped.
      #   start$arcs           arcs already in the graph.
      # !is.listed(whitelist)  exclude whitelisted arcs.
      if (!is.null(blacklist))
        to.be.dropped = start$arcs[!is.row.equal(start$arcs, whitelist), , drop = FALSE]
      else
        to.be.dropped = start$arcs

      # if there is any arc in the graph, try to remove it.
      if (nrow(to.be.dropped) > 0) {

        if (debug) {

          cat("----------------------------------------------------------------\n")
          cat("* trying to remove one of", nrow(start$arcs), "arcs.\n")

        }#THEN

        apply(to.be.dropped, 1, drop.step)

      }#THEN
      else if (debug) {

        cat("----------------------------------------------------------------\n")
        cat("* no arc to be removed.\n")

      }#THEN

      # create a list of arcs to be reversed.
      if (!is.null(blacklist))
        to.be.reversed = start$arcs[!is.row.equal(start$arcs[, c(2,1), drop = FALSE], blacklist), , drop = FALSE]
      else
        to.be.reversed = start$arcs

      if (nrow(to.be.reversed) > 0) {

        if (debug) {

          cat("----------------------------------------------------------------\n")
          cat("* trying to reverse one of", nrow(start$arcs), "arcs.\n")

        }#THEN

        apply(to.be.reversed, 1, reverse.step)

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
    if (identical(start, end)) {

      # the likelihood score has this nasty habit to be too small for its own
      # good, and is often rounded to zero due to machine precision limit. If
      # that's the case, switch to the loglikelihood.
      if ((score == "lik") && (sum(reference.score) == 0)) {

        warning("switching to loglikelihood due to machine precision limits.")

        score = "loglik"
        reference.score = sapply(names(start$nodes), loglik.node,
                            x = start, data = x)

      }#THEN
      else
        break

    }#THEN

    # update the cached values of the end network.
    end$nodes = cache.structure(nodes, end$arcs)
    # reset the score delta.
    end$score.delta = 0
    # update the starting network for the next iteration.
    start = end
    # reset the reference score.
    if (!is.null(end$updates))
      reference.score[names(end$updates)] = end$updates

    if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("* current network is :\n")
      print(start)
      cat("* current score:", sum(reference.score), "\n")

    }#THEN

  }#REPEAT

  # remove score.delta and updated elemnts from the return value.
  end$updates = end$score.delta = NULL

  end

}#HILL.CLIMBING

hill.climbing.optimized = function(x, start, whitelist, blacklist, score,
    extra.args, debug) {

  # cache nodes' labels.
  nodes = colnames(x)
  # cache the number of nodes.
  n.nodes = length(nodes)
  # create the prospective model.
  end = start
  # set the score delta.
  end$score.delta = 0
  # check whether the score is score-equivalent.
  score.equivalence = score %in% score.equivalent.scores

  if (score.equivalence) {

    # set up the score-equivalence 'shadow' cache.
    shadow.cache.add = data.frame(from = character(0), to = character(0),
                         delta = numeric(0), stringsAsFactors = FALSE)

  }#THEN

  # function to add an arc to the graph.
  add.step = function(arc) {

    # initialize the buffered score delta to NA.
    w = NA

    # there's no need to call has.path() here; in most cases the arc is not
    # going to be added at all (so it would be pointless).

    # the cache is empty if the network has not yet been updated; look for
    # a score equivalent arc whose score delta is already cached.
    if (is.null(start$updates) && score.equivalence) {

       # retrieve the last "known good" score delta from this other cache.
        w = shadow.cache.add[(shadow.cache.add[,1] == arc[1]) & (shadow.cache.add[,2] == arc[2]), "delta"]

        if (!is.na(w) && (length(w) > 0)) {

          if (debug) cat("  > found", arc[2], "->", arc[1], "score equivalent to",
                       arc[1], "->", arc[2], "in the shadow cache.\n")

        }#THEN

    }#THEN

    # if there is a valid score cache, look for a cached score.
    if (!is.null(start$updates)) {

      # the node with the incoming arrow must have been updated in
      # the last iteration; otherwise is cached.
      if (!(arc[2] %in% names(start$updates))) {

        # retrieve the last "known good" score delta for this arc.
        w = score.cache.add[(score.cache.add[,1] == arc[1]) & (score.cache.add[,2] == arc[2]), "delta"]

      }#THEN

    }#THEN

    # if the score delta is there ...
    if (!is.na(w) && (length(w) > 0)) {

      # ... and is negative do nothing.
      if (w <= 0) {

        if (debug) {

          cat("  > trying to add", arc[1], "->", arc[2], "(cached!).\n")
          cat("    > delta between scores for nodes", arc, "is", w, ".\n")

        }#THEN

        return(w)

      }#THEN
      else if ((w > 0) && (w > end$score.delta)) {

        # if the arc creates cycles, do not add it.
        if (has.path(arc[2], arc[1], nodes, amat)) {

          if (debug) cat("  > trying to add", arc[1], "->", arc[2], "(cycles!).\n")

          return(NA)

        }#THEN

        if (debug) {

          cat("  > trying to add", arc[1], "->", arc[2], "(cached!).\n")
          cat("    > delta between scores for nodes", arc, "is", w, ".\n")
          cat("    @ adding", arc[1], "->", arc[2], ".\n")

        }#THEN

        # no need to call set.arc.direction, the arc is not present.
        start$arcs = rbind(start$arcs, arc, deparse.level = 0)
        start$score.delta = w
        start$updates = reference.score[arc[2]] + w
        names(start$updates) = arc[2]

        assign("end", start, envir = sys.frame(-3))

        return(w)

      }#THEN

    }#THEN

    # if the arc creates cycles, do not add it.
    if (has.path(arc[2], arc[1], nodes, amat)) {

      if (debug) cat("  > trying to add", arc[1], "->", arc[2], "(cycles!).\n")

      return(NA)

    }#THEN

    if (debug)
      cat("  > trying to add", arc[1], "->", arc[2], ".\n")

    # compare the scores of the two networks.
    better = score.delta(arc = arc, network = start, data = x,
               score = score, score.delta = end$score.delta,
               reference.score = reference.score, op = "set",
               extra = extra.args, debug = debug)

    if (better$bool) {

      if (debug) cat("    @ adding", arc[1], "->", arc[2], ".\n")

      # no need to call set.arc.direction, the arc is not present.
      start$arcs = rbind(start$arcs, arc, deparse.level = 0)
      start$score.delta = better$delta
      start$updates = better$updates

      # if score.cache.add is used, thie apply() is called inside data.frame(),
      # so the right scope is one step further than usual.
      assign("end", start, envir = sys.frame(-3))

    }#THEN

   # update the score-equivalence 'shadow' cache until I have a real cache.
   if (score.equivalence && is.null(start$updates)) {

     assign("shadow.cache.add", rbind(shadow.cache.add, data.frame(from = arc[2],
       to = arc[1], delta = better$delta, stringsAsFactors = FALSE)),
       envir = sys.frame(-3))

   }#THEN

    return(better$delta)

  }#ADD.STEP

  # function to drop an arc from the graph.
  drop.step = function(arc) {

    # the cache is empty if the network has not yet been updated.
    if (!is.null(start$updates)) {

      # the node with the incoming arrow must have been updated in
      # the last iteration; otherwise is cached.
      if (!(arc[2] %in% names(start$updates))) {

        # retrieve the last "known good" score delta for this arc.
        w = score.cache.drop[(score.cache.drop[,1] == arc[1]) & (score.cache.drop[,2] == arc[2]), "delta"]

        # if the score delta is there ...
        if (!is.na(w) && (length(w) > 0)) {

          # ... and is negative do nothing.
          if (w <= 0) {

            if (debug) {

              cat("  > trying to remove", arc[1], "->", arc[2], "(cached!).\n")
              cat("    > delta between scores for nodes", arc, "is", w, ".\n")

            }#THEN

            return(w)

          }#THEN
          else if ((w > 0) && (w > end$score.delta)) {

            if (debug) {

              cat("  > trying to remove", arc[1], "->", arc[2], "(cached!).\n")
              cat("    > delta between scores for nodes", arc, "is", w, ".\n")
              cat("    @ removing", arc[1], "->", arc[2], ".\n")

            }#THEN

            # drop the arc from the network.
            start$arcs = drop.arc.backend(start$arcs, arc)
            start$score.delta = w
            start$updates = reference.score[arc[2]] + w
            names(start$updates) = arc[2]

            assign("end", start, envir = sys.frame(-3))

            return(w)

          }#THEN

        }#THEN


      }#THEN

    }#THEN

    if (debug)
      cat("  > trying to remove", arc[1], "->", arc[2], ".\n")

    # compare the scores of the two networks.
    better = score.delta(arc = arc, network = start, data = x,
               score = score, score.delta = end$score.delta,
               reference.score = reference.score, op = "drop",
               extra = extra.args, debug = debug)

    if (better$bool) {

      if (debug) cat("    @ removing", arc[1], "->", arc[2], ".\n")

      # drop the arc from the network.
      start$arcs = drop.arc.backend(start$arcs, arc)
      start$score.delta = better$delta
      start$updates = better$updates

      # if score.cache.drop is used, thie apply() is called inside
      # data.frame(), so the right scope is one step further than usual.
      assign("end", start, envir = sys.frame(-3))

    }#THEN

    return(better$delta)

  }#DROP.STEP

  # function to reverse an arc in the graph.
  reverse.step = function(arc) {

    # the cache is empty if the network has not yet been updated.
    if (!is.null(start$updates)) {

      # the node with the incoming arrow must have been updated in
      # the last iteration; otherwise is cached.
      if (!any(arc %in% names(start$updates))) {

        # retrieve the last "known good" score delta for this arc.
        w = score.cache.reverse[(score.cache.reverse[,1] == arc[1]) & (score.cache.reverse[,2] == arc[2]), "delta"]

        # if the score delta is there ...
        if (!is.na(w) && (length(w) > 0)) {

          # ... and is negative do nothing.
          if (w <= 0) {

            if (debug) {

              cat("  > trying to reverse", arc[1], "->", arc[2], "(cached!).\n")
              cat("    > delta between scores for nodes", arc, "is", w, ".\n")

            }#THEN

            return(w)

          }#THEN

        }#THEN

      }#THEN

    }#THEN

    # if the arc creates cycles, do not reverse it.
    if (has.path(arc[1], arc[2], nodes, amat, exclude.direct = TRUE)) {

      if (debug) cat("  > trying to reverse", arc[1], "->", arc[2], "(cycles!).\n")

      return(FALSE)

    }#THEN

    if (debug)
      cat("  > trying to reverse", arc[1], "->", arc[2], ".\n")

    # compare the scores of the two networks.
    better = score.delta(arc = arc, network = start, data = x,
               score = score, score.delta = end$score.delta,
               reference.score = reference.score, op = "reverse",
               extra = extra.args, debug = debug)

    if (better$bool) {

      if (debug) cat("    @ reversing", arc[1], "->", arc[2], ".\n")

      start$arcs = reverse.arc.backend(arc[1], arc[2], start$arcs)
      start$score.delta = better$delta
      start$updates = better$updates

      # if score.cache.reverse is used, thie apply() is called inside
      # data.frame(), so the right scope is one step further than usual.
      assign("end", start, envir = sys.frame(-3))

    }#THEN

    return(better$delta)

  }#REVERSE.STEP

  # set the reference score.
  if (score == "k2") {

    reference.score = sapply(names(start$nodes), dirichlet.node,
                        x = start, data = x)

  }#THEN
  else if (score %in% c("bde", "dir")) {

    reference.score = sapply(names(start$nodes), dirichlet.node,
                        x = start, imaginary.sample.size = extra.args$iss,
                        data = x)

  }#THEN
  else if (score %in% c("lik", "loglik")) {

    reference.score = sapply(names(start$nodes), loglik.node,
                        x = start, data = x)
    if (score == "lik") reference.score = exp(reference.score)

  }#THEN
  else if (score %in% c("aic", "bic")) {

    reference.score = sapply(names(start$nodes), aic.node, x = start, data = x,
                        k = extra.args$k)

  }#THEN
  else if (score == "bge") {

    reference.score = sapply(names(start$nodes), bge.node,
                        x = start, imaginary.sample.size = extra.args$iss,
                        phi = extra.args$phi, data = x)

  }#THEN

  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* starting from the following network:\n")
    print(start)
    cat("* current score:", sum(reference.score), "\n")
    cat("* whitelisted arcs are:\n")
    if (!is.null(whitelist)) print(whitelist)
    cat("* blacklisted arcs are:\n")
    if (!is.null(blacklist)) print(blacklist)

  }#THEN

  if (!is.null(blacklist))
    blmat = arcs2amat(blacklist, nodes)

  repeat {

    amat = arcs2amat(start$arcs, nodes)

    # create the list of arcs to be added.
    #   1               all the possibile arcs.
    # - amat            exclude arcs already in the graph.
    # - t(amat)         exclude the reverse of those arcs.
    # - diag(n.nodes)   exclude self-loops.
    # - blmat		exclude blacklisted arcs.
    if (!is.null(blacklist))
      to.be.added = amat2arcs(1 - amat - t(amat) - diag(n.nodes) - blmat, nodes)
    else
      to.be.added = amat2arcs(1 - amat - t(amat) - diag(n.nodes), nodes)

    if (nrow(to.be.added) > 0) {

      if (debug) {

        cat("----------------------------------------------------------------\n")
        cat("* trying to add one of", nrow(to.be.added), "arcs.\n")

      }#THEN

      # try to add any available arc.
      score.cache.add = data.frame(to.be.added,
          delta = apply(to.be.added, 1, add.step), stringsAsFactors = FALSE)

    }#THEN
    else if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("* no arc to be added.\n")

    }#THEN

    if (nrow(start$arcs) > 0) {

      # create a list of arcs to be dropped.
      #   start$arcs           arcs already in the graph.
      # !is.listed(whitelist)  exclude whitelisted arcs.
      if (!is.null(blacklist))
        to.be.dropped = start$arcs[!is.row.equal(start$arcs, whitelist), , drop = FALSE]
      else
        to.be.dropped = start$arcs

      # if there is any arc in the graph, try to remove it.
      if (nrow(to.be.dropped) > 0) {

        if (debug) {

          cat("----------------------------------------------------------------\n")
          cat("* trying to remove one of", nrow(start$arcs), "arcs.\n")

        }#THEN

        # try to remove any arc in the graph.
        score.cache.drop = data.frame(to.be.dropped,
            delta = apply(to.be.dropped, 1, drop.step), stringsAsFactors = FALSE)

      }#THEN
      else if (debug) {

        cat("----------------------------------------------------------------\n")
        cat("* no arc to be removed.\n")

      }#THEN

      # create a list of arcs to be reversed.
      if (!is.null(blacklist))
        to.be.reversed = start$arcs[!is.row.equal(start$arcs[, c(2,1), drop = FALSE], blacklist), , drop = FALSE]
      else
        to.be.reversed = start$arcs

      if (nrow(to.be.reversed) > 0) {

        if (debug) {

          cat("----------------------------------------------------------------\n")
          cat("* trying to reverse one of", nrow(start$arcs), "arcs.\n")

        }#THEN

        # try to remove any arc in the graph.
        score.cache.reverse = data.frame(to.be.reversed,
            delta = apply(to.be.reversed, 1, reverse.step), stringsAsFactors = FALSE)

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
    if (identical(start, end)) {

      # the likelihood score has this nasty habit to be too small for its own
      # good, and is often rounded to zero due to machine precision limit. If
      # that's the case, switch to the loglikelihood.
      if ((score == "lik") && (sum(reference.score) == 0)) {

        warning("switching to loglikelihood due to machine precision limits.")

        score = "loglik"
        reference.score = sapply(names(start$nodes), loglik.node,
                            x = start, data = x)

      }#THEN
      else
        break

    }#THEN

    # update the cached values of the end network.
    end$nodes = cache.structure(nodes, end$arcs)
    # reset the score delta.
    end$score.delta = 0
    # update the starting network for the next iteration.
    start = end
    # reset the reference score.
    if (!is.null(end$updates))
      reference.score[names(end$updates)] = end$updates

    if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("* current network is :\n")
      print(start)
      cat("* current score:", sum(reference.score), "\n")

    }#THEN

  }#REPEAT

  # remove score.delta and updated elemnts from the return value.
  end$updates = end$score.delta = NULL

  end

}#HILL.CLIMBING

