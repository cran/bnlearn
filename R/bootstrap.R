
# simple parametric and nonparametric bootstrap implementation.
bootstrap.backend = function(data, statistic, R, m, sim = "ordinary",
    algorithm, algorithm.args = list(), statistic.args = list(),
    debug = FALSE) {

  # allocate the result list.
  res = vector(R, mode = "list")

  # initialize the bayesian network used by the paramentric bootstrap.
  if (sim == "parametric") {

    if (algorithm %in% always.dag.result) {

      net = do.call(algorithm, c(list(x = data), algorithm.args))

      if (debug) {

        cat("* initial network for parametric bootstrap is:\n")
        print(net)

      }#THEN

    }#THEN
    else {

      stop(paste("this learning algorithm may result in a partially",
             "directed dag, which is not handled by parametric bootstrap."))

    }#ELSE

  }#THEN

  for (r in seq_len(R)) {

    if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("* bootstrap replicate", r, ".\n")

    }#THEN

    # generate the r-th bootstrap sample.
    if (sim == "ordinary")
      replicate = data[sample(nrow(data), m, replace = TRUE), , drop = FALSE]
    else if (sim == "parametric")
      if (is.data.discrete(data))
        replicate = rbn.discrete(x = net, n = m, data = data)
      else
        replicate = rbn.continuous(x = net, n = m, data = data)

    if (debug)
      cat("* learning bayesian network structure.\n")

    # learn the network structure from the bootstrap sample.
    net = do.call(algorithm, c(list(x = replicate), algorithm.args))

    if (debug) {

      print(net)
      cat("* applying user-defined statistic.\n")

    }#THEN

    # apply the user-defined function to the newly-learned bayesian network;
    # the bayesian network is passed as the first argument hoping it will end
    # at the right place thanks to the positional matching.
    res[[r]] = do.call(statistic, c(list(net), statistic.args))

    if (debug) {

      cat("  > the function returned:\n")
      print(res[[r]]);

    }#THEN

  }#FOR

  return(res)

}#BOOTSTRAP.BACKEND

