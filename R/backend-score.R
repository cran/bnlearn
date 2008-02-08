
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

      return(FALSE)

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

# compare two network scores in an efficient way.
score.delta = function(arc, network, data, score, score.delta, reference.score, op, extra, debug) {

  if (op == "reverse") {

    # do a minimal update to the network structure.
    network$nodes[[arc[1]]]$parents = c(network$nodes[[arc[1]]]$parents, arc[2])
    network$nodes[[arc[2]]]$parents =
      network$nodes[[arc[2]]]$parents[network$nodes[[arc[2]]]$parents != arc[1]]

    if (score == "k2") {

      new.score = c(dirichlet.node(node = arc[1], x = network, data = data),
                  dirichlet.node(node = arc[2], x = network, data = data))

    }#THEN
    else if (score %in% c("bde", "dir")) {

      new.score = c(dirichlet.node(node = arc[1], x = network, imaginary.sample.size = extra$iss, data = data),
                  dirichlet.node(node = arc[2], x = network, imaginary.sample.size = extra$iss, data = data))

    }#THEN
    if (score %in% c("lik", "loglik")) {

      new.score = c(loglik.node(node = arc[1], x = network, data = data),
                  loglik.node(node = arc[2], x = network, data = data))
      if (score == "lik") new.score = exp(new.score)

    }#THEN
    else if (score %in% c("aic", "bic")) {

      new.score = c(aic.node(node = arc[1], x = network, k = extra$k, data = data),
                  aic.node(node = arc[2], x = network, k = extra$k, data = data))

    }#THEN
    else if (score == "bge") {

      new.score = c(bge.node(node = arc[1], x = network, imaginary.sample.size = extra$iss, phi = extra$phi, data = data),
                  bge.node(node = arc[2], x = network, imaginary.sample.size = extra$iss, phi = extra$phi, data = data))

    }#THEN

    # set the names of the updated elements.
    names(new.score) = arc

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

    # compute the new score for the node arc[2].
    if (score == "k2") {

      new.score = dirichlet.node(node = arc[2], x = network, data = data)

    }#THEN
    else if (score %in% c("bde", "dir")) {

      new.score = dirichlet.node(node = arc[2], x = network, imaginary.sample.size = extra$iss, data = data)

    }#THEN
    if (score %in% c("lik", "loglik")) {

      new.score = loglik.node(node = arc[2], x = network, data = data)
      if (score == "lik") new.score = exp(new.score)

    }#THEN
    else if (score %in% c("aic", "bic")) {

      new.score = aic.node(node = arc[2], x = network, k = extra$k, data = data)

    }#THEN
    else if (score == "bge") {

      new.score = bge.node(node = arc[2], x = network, imaginary.sample.size = extra$iss, phi = extra$phi, data = data)

    }#THEN

    # set the names of the updated elements.
    names(new.score) = arc[2]

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

