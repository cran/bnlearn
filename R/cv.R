
crossvalidation = function(data, bn, loss = NULL, k = 5, algorithm.args,
    loss.args, fit, fit.args, cluster = NULL, debug = FALSE) {

  n = nrow(data)

  # shuffle the data to get unbiased splits.
  kcv = split(sample(n), seq_len(k))
  # store the length of each test set.
  kcv.length = sapply(kcv, length)

  if (debug) {

    cat("* splitting", n, "data in", k, "subsets.\n")
    cat("----------------------------------------------------------------\n")

  }#THEN

  if (is.character(bn)) {

    if (!is.null(cluster)) {

      kcv = parLapply(cluster, kcv, bn.cv.algorithm, data = data,
              algorithm = bn, algorithm.args = algorithm.args, loss = loss,
              loss.args = loss.args, fit = fit, fit.args = fit.args,
              debug = debug)

    }#THEN
    else {

      kcv = lapply(kcv, bn.cv.algorithm, data = data, algorithm = bn,
              algorithm.args = algorithm.args, loss = loss,
              loss.args = loss.args, fit = fit, fit.args = fit.args,
              debug = debug)

    }#ELSE

  }#THEN
  else if (is(bn, "bn")) {

    if (!is.null(cluster)) {

      kcv = parLapply(cluster, kcv, bn.cv.structure, data = data, bn = bn,
              loss = loss, loss.args = loss.args, fit = fit, fit.args = fit.args,
              debug = debug)
    }#THEN
    else {

      kcv = lapply(kcv, bn.cv.structure, data = data, bn = bn, loss = loss,
              loss.args = loss.args, fit = fit, fit.args = fit.args,
              debug = debug)

    }#ELSE

  }#THEN

  if (debug) {

    cat("* summary of the observed values for the loss function:\n")
    print(summary(unlist(sapply(kcv, '[', 'loss'))))

  }#THEN

  # compute the mean of the observed values of the loss function, weighted
  # to account for unequal-length splits.
  mean = weighted.mean(unlist(sapply(kcv, '[', 'loss')), kcv.length)

  # reset the names of the elements of the return value.
  names(kcv) = NULL
  # add some useful attributes to the return value.
  kcv = structure(kcv, class = "bn.kcv", loss = loss, args = loss.args,
          bn = bn, mean = mean)

  return(kcv)

}#BN.CV

bn.cv.algorithm = function(test, data, algorithm, algorithm.args, loss,
    loss.args, fit, fit.args, debug = FALSE) {

  if (debug)
    cat("* learning the structure of the network from the training sample.\n")

  # learning the structure of the network.
  net = do.call(algorithm, c(list(x = data[-test, ]), algorithm.args))

  if (debug)
    print(net)

  # go on with fitting the parameters.
  bn.cv.structure(test = test, data = data, bn = net, loss = loss,
    loss.args = loss.args, fit = fit, fit.args = fit.args, debug = debug)

}#BN.CV.ALGORITHM

bn.cv.structure = function(test, data, bn, loss, loss.args, fit, fit.args,
    debug = FALSE) {

  if (debug)
    cat("* fitting the parameters of the network from the training sample.\n")

  # use score equivalence to compute log-likelihood losses when the network
  # returned by the learning algorithm is a CPDAG; extend it to a DAG (which
  # has the same log-likelihoos because it's in the same equivalence class)
  # and use the result in place of the original network.
  if ((loss %in% c("logl", "logl-g")) &&
      !is.dag(arcs = bn$arcs, nodes = names(bn$nodes))) {

    bn = cpdag.extension(cpdag.backend(bn))

  }#THEN

  # check the extra arguments.
  fit.args = check.fitting.args(fit, bn, data[-test, ], fit.args)
  # fit the parameters.
  net = bn.fit.backend(x = bn, data = data[-test, ], method = fit,
          extra.args = fit.args)

  if (debug)
    cat("* applying the loss function to the data from the test sample.\n")

  obs.loss = loss.function(fitted = net, data = data[test, ], loss = loss,
               extra.args = loss.args, debug = debug)

  if (debug)
    cat("----------------------------------------------------------------\n")

  return(c(list(test = test, fitted = net), obs.loss))

}#BN.CV.STRUCTURE

