
# generic frontend to {non,}parametric bootstrap.
bn.boot = function(data, statistic, R = 200, m = nrow(data), algorithm,
    algorithm.args = list(), statistic.args = list(), cluster = NULL,
    debug = FALSE) {

  # check the data are there.
  check.data(data)
  # check the number of bootstrap replicates.
  R = check.replicates(R)
  # check the size of each bootstrap sample.
  m = check.bootsize(m, data)
  # check debug.
  check.logical(debug)
  # check the learning algorithm.
  check.learning.algorithm(algorithm)
  # check the extra arguments for the learning algorithm.
  algorithm.args = check.learning.algorithm.args(algorithm.args)
  # check the custom statistic function.
  statistic = match.fun(statistic)
  # if the statistic fuction is I(), replace it with a NOP to avoid troubles
  # with class dispatch.
  if (isTRUE(all.equal(statistic, base::I)))
    statistic = function(x) x

  # check the extra arguments for the statistic function.
  if (!is.list(statistic.args))
    statistic.args = as.list(statistic.args)

  # check the cluster.
  if (!is.null(cluster)) {

    check.cluster(cluster)

    # disable debugging, the slaves do not cat() here.
    if (debug) {

      warning("disabling debugging output for parallel computing.")
      debug = FALSE

    }#THEN

  }#THEN

  bootstrap.backend(data = data, statistic = statistic, R = R, m = m,
    algorithm = algorithm, algorithm.args = algorithm.args,
    statistic.args = statistic.args, cluster = cluster, debug = debug)

}#BNBOOT

# compute arcs' strength via nonparametric bootstrap.
boot.strength = function(data, cluster = NULL, R = 200, m = nrow(data),
    algorithm, algorithm.args = list(), cpdag = TRUE, debug = FALSE) {

  # check the data are there.
  data.info = check.data(data)
  # check the number of bootstrap replicates.
  R = check.replicates(R)
  # check the size of each bootstrap sample.
  m = check.bootsize(m, data)
  # check debug.
  check.logical(debug)
  # check cpdag.
  check.logical(cpdag)
  # check the learning algorithm.
  check.learning.algorithm(algorithm)
  # check the extra arguments for the learning algorithm.
  algorithm.args = check.learning.algorithm.args(algorithm.args)

  # check the cluster.
  if (!is.null(cluster)) {

    check.cluster(cluster)

    # set up the slave processes.
    slaves.setup(cluster)
    # disable debugging, the slaves do not cat() here.
    if (debug) {

      warning("disabling debugging output for parallel computing.")
      debug = FALSE

    }#THEN

  }#THEN

  res = arc.strength.boot(data = data, cluster = cluster, R = R, m = m,
          algorithm = algorithm, algorithm.args = algorithm.args, arcs = NULL,
          cpdag = cpdag, debug = debug)

  # add extra information for strength.plot().
  res = structure(res, method = "bootstrap", threshold = threshold(res),
          class = c("bn.strength", class(res)))
  if (data.info$type == "mixed-cg")
    attr(res, "illegal") = list.cg.illegal.arcs(names(data), data)

  return(res)

}#BOOT.STRENGTH

# perform cross-validation.
bn.cv = function(data, bn, loss = NULL, ...,
    algorithm.args = list(), loss.args = list(), fit = "mle",
    fit.args = list(), method = "k-fold", cluster = NULL, debug = FALSE) {

  # check the data are there.
  data.info = check.data(data)
  nodes = names(data)
  # check the cross-validation method.
  method = check.cv.method(method)
  # check the extra arguments for the cross-validation method.
  extra.args = check.cv.args(method, list(...), data)
  # check the fitting method (maximum likelihood, bayesian, etc.)
  check.fitting.method(fit, data)
  # check the cluster.
  if (!is.null(cluster)) {

    check.cluster(cluster)

    # disable debugging, the slaves do not cat() here.
    if (debug) {

      warning("disabling debugging output for parallel computing.")
      debug = FALSE

    }#THEN

  }#THEN

  if (is.character(bn)) {

    # check the learning algorithm.
    check.learning.algorithm(bn)
    # check the loss function.
    loss = check.loss(loss, data, bn)
    # check the extra arguments for the learning algorithm.
    algorithm.args = check.learning.algorithm.args(algorithm.args)
    # since we have no other way to guess, copy the label of the target
    # variable from the parameters of the classifier.
    if ((loss %in% c("pred", "pred-exact", "pred-lw")) &&
        (is.null(loss.args$target)) && (bn %in% classifiers))
      loss.args$target = algorithm.args$training
    # check the extra arguments passed down to the loss function.
    loss.args = check.loss.args(loss, bn, nodes, data, loss.args)

  }#THEN
  else if (is(bn, "bn")) {

    if (!identical(algorithm.args, list()))
      warning("no learning algorithm is used, so 'algorithm.args' will be ignored.")
    # check whether the data agree with the bayesian network.
    check.bn.vs.data(bn, data)
    # check the loss function.
    loss = check.loss(loss, data, bn)
    # it is impossible to learn the parameters if the network structure is only
    # partially directed and cannot be extended to a completely directed graph.
    if (!is.dag(arcs = bn$arcs, nodes = nodes)) {

      # trying to extend a skeleton (instead of a CPDAG) is probably not
      # meaningful.
      if (!is.null(bn$learning$undirected) && bn$learning$undirected)
        warning(deparse(substitute(bn)), " is just a skeleton (no arc directions ",
          "have been learned) and trying to extend it is probably wrong.")

      # try to extend the network into a completely directed graph, keeping as
      # many of the original arc directions are possible.
      bn = cpdag.extension(cpdag.backend(bn, moral = TRUE, wlbl = TRUE))

      if (!is.dag(arcs = bn$arcs, nodes = nodes))
        stop("the network from the 'bn' argument has no consistent extension.")

    }#THEN

    # check bn.naive objects if any.
    if (is(bn, "bn.naive"))
      check.bn.naive(bn)
    # check the extra arguments passed down to the loss function.
    loss.args = check.loss.args(loss, bn, nodes, data, loss.args)

  }#THEN
  else {

    stop("'bn' must be either the label of a learning algorithm or a bn object.")

  }#ELSE

  # allocate and populate the return value.
  actual.runs = ifelse(is.null(extra.args$runs), 1, extra.args$runs)

  result = structure(vector(actual.runs, mode = "list"), class = "bn.kcv.list")

  for (r in seq(actual.runs))
    result[[r]] = crossvalidation(data = data, bn = bn, loss = loss,
                    k = extra.args$k, m = extra.args$m,
                    folds = extra.args$folds[[r]],
                    algorithm.args = algorithm.args, loss.args = loss.args,
                    fit = fit, fit.args = fit.args, method = method,
                    cluster = cluster, data.info = data.info, debug = debug)

  # return a bn.kcv object (for a single run) or a bn.kcv.list object (for
  # multiple runs).
  if (actual.runs == 1)
    return(result[[1]])
  else
    return(result)

}#BN.CV

# extract loss values from bn.kcv and bn.kcv.list objects.
loss = function(x) {

  if (is(x, "bn.kcv"))
    losses = attr(x, "mean")
  else if (is(x, "bn.kcv.list"))
    losses = sapply(x, function(x) attr(x, "mean"))
  else
    stop("x must be an object of class 'bn.kcv' or 'bn.kcv.list'.")

  return(losses)

}#LOSS
