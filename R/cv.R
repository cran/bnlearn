
crossvalidation = function(data, bn, loss = NULL, k = 10,
    m = ceiling(nrow(data)/10), folds, algorithm.args, loss.args, fit,
    fit.args, method, cluster = NULL, data.info, debug = FALSE) {

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
  else if (method == "custom-folds") {

    # the folds are the custom folds specified by the user.
    kcv = folds
    # store the length of each test set.
    kcv.length = sapply(kcv, length)

  }#THEN

  if (is.character(bn)) {

    if (!is.null(cluster)) {

      kcv = parallel::parLapplyLB(cluster, kcv, bn.cv.algorithm, data = data,
              algorithm = bn, algorithm.args = algorithm.args, loss = loss,
              loss.args = loss.args, fit = fit, fit.args = fit.args,
              data.info = data.info, debug = debug)

    }#THEN
    else {

      kcv = lapply(kcv, bn.cv.algorithm, data = data, algorithm = bn,
              algorithm.args = algorithm.args, loss = loss,
              loss.args = loss.args, fit = fit, fit.args = fit.args,
              data.info = data.info, debug = debug)

    }#ELSE

  }#THEN
  else if (is(bn, "bn")) {

    if (!is.null(cluster)) {

      kcv = parallel::parLapplyLB(cluster, kcv, bn.cv.structure, data = data,
              bn = bn, loss = loss, loss.args = loss.args, fit = fit,
              fit.args = fit.args, data.info = data.info, debug = debug)

    }#THEN
    else {

      kcv = lapply(kcv, bn.cv.structure, data = data, net = bn, loss = loss,
              loss.args = loss.args, fit = fit, fit.args = fit.args,
              data.info = data.info, debug = debug)

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
    loss.args, fit, fit.args, data.info, debug = FALSE) {

  if (debug)
    cat("* learning the structure of the network from the training sample.\n")

  # learning the structure of the network.
  net = do.call(algorithm, c(list(x = data[-test, ]), algorithm.args))

  if (debug)
    print(net)

  # go on with fitting the parameters.
  bn.cv.structure(test = test, data = data, net = net, loss = loss,
    loss.args = loss.args, fit = fit, fit.args = fit.args,
    data.info = data.info, debug = debug)

}#BN.CV.ALGORITHM

bn.cv.structure = function(test, data, net, loss, loss.args, fit, fit.args,
    data.info, debug = FALSE) {

  if (debug)
    cat("* fitting the parameters of the network from the training sample.\n")

  # if the network is not completely directed, try to extend it before moving on
  # and performing parameter learning.
  if (!is.dag(arcs = net$arcs, nodes = names(net$nodes))) {

    # trying to extend a skeleton (instead of a CPDAG) is probably not
    # meaningful.
    if (!is.null(net$learning$undirected) && net$learning$undirected)
      warning("the network in one of the folds is just a skeleton (no arc ",
        "directions have been learned) and trying to extend it is probably wrong.")

    net = cpdag.extension(cpdag.backend(net, moral = TRUE, wlbl = TRUE))

    if (loss %in% c("pred", "cor", "mse")) {

      # for some loss functions just identifying the parents of the target node
      # is enough.
      target = loss.args$target

      if (!setequal(net$nodes[[target]]$nbr,
            union(net$nodes[[target]]$parents, net$nodes[[target]]$children)))
        stop("undirected arcs around the target node ", target, ".")

    }#THEN
    else {

      # other loss functions require the network to be completely directed, so
      # so that all local distributions can be estimated.
      if (any(which.undirected(net$arcs, names(net$nodes))))
        stop("no consistent extension for the network in one of the folds.")

    }#ELSE

  }#THEN

  # check the extra arguments.
  fit.args = check.fitting.args(fit, net, data[-test, ], fit.args)
  # fit the parameters.
  fitted = bn.fit.backend(x = net, data = data[-test, ], method = fit,
             extra.args = fit.args, data.info = data.info)

  # in the case of naive Bayes and TAN models, the prior must be computed on
  # the training sample for each fold to match the behaviour of the default
  # for non-cross-validated models.
  if (is(fitted, c("bn.naive", "bn.tan")))
    if (all(loss.args$prior == 1))
      loss.args$prior = fitted[[attr(fitted, "training")]]$prob

  if (debug)
    cat("* applying the loss function to the data from the test sample.\n")

  obs.loss = loss.function(fitted = fitted, data = data[test, ], loss = loss,
               extra.args = loss.args, debug = debug)

  if (debug)
    cat("----------------------------------------------------------------\n")

  return(c(list(test = test, fitted = fitted, learning = net$learning), obs.loss))

}#BN.CV.STRUCTURE

