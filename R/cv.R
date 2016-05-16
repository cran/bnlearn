
crossvalidation = function(data, bn, loss = NULL, k = 10,
    m = ceiling(nrow(data)/10), algorithm.args, loss.args, fit, fit.args,
    method, cluster = NULL, debug = FALSE) {

  n = nrow(data)

  if (method == "k-fold") {

    # shuffle the data to get unbiased splits (do not warn about fold size).
    kcv = suppressWarnings(split(sample(n), seq_len(k)))
    # store the length of each test set.
    kcv.length = sapply(kcv, length)

    if (debug) {

      cat("* splitting", n, "data in", k, "subsets.\n")
      cat("----------------------------------------------------------------\n")

    }#THEN

  }#THEN
  else if (method == "hold-out") {

    # sample m observations without replacement.
    kcv = lapply(seq(k), function(x) sample(n, m))
    # all test sets have the same length, a dummy works just fine.
    kcv.length = rep(m, k)

  }#THEN

  if (is.character(bn)) {

    if (!is.null(cluster)) {

      kcv = parallel::parLapply(cluster, kcv, bn.cv.algorithm, data = data,
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

      kcv = parallel::parLapply(cluster, kcv, bn.cv.structure, data = data,
              bn = bn, loss = loss, loss.args = loss.args, fit = fit,
              fit.args = fit.args, debug = debug)

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

  # aggregate the loss functions computed over the k folds.
  mean = kfold.loss.postprocess(kcv, kcv.length, loss, loss.args, data)

  # reset the names of the elements of the return value.
  names(kcv) = NULL
  # add some useful attributes to the return value.
  kcv = structure(kcv, class = "bn.kcv", loss = loss, args = loss.args,
          bn = bn, mean = mean, method = method)

  return(kcv)

}#CROSSVALIDATION

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
  if (loss %in% c("logl", "logl-g", "logl-cg") &&
      !is.dag(arcs = bn$arcs, nodes = names(bn$nodes))) {

    bn = cpdag.extension(cpdag.backend(bn, wlbl = TRUE))

  }#THEN

  # check the extra arguments.
  fit.args = check.fitting.args(fit, bn, data[-test, ], fit.args)
  # fit the parameters.
  net = bn.fit.backend(x = bn, data = data[-test, ], method = fit,
          extra.args = fit.args)

  # in the case of naive Bayes and TAN models, the prior must be computed on
  # the training sample for each fold to match the behaviour of the default
  # for non-cross-validated models.
  if (is(net, c("bn.naive", "bn.tan")))
    if (all(loss.args$prior == 1))
      loss.args$prior = net[[attr(net, "training")]]$prob

  if (debug)
    cat("* applying the loss function to the data from the test sample.\n")

  obs.loss = loss.function(fitted = net, data = data[test, ], loss = loss,
               extra.args = loss.args, debug = debug)

  if (debug)
    cat("----------------------------------------------------------------\n")

  return(c(list(test = test, fitted = net), obs.loss))

}#BN.CV.STRUCTURE

