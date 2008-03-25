
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

