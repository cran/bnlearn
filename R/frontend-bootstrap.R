
# generic frontend to {non,}parametric bootstrap.
bn.boot = function(data, statistic, R = 200, m = nrow(data), sim = "ordinary",
    algorithm, algorithm.args = list(), statistic.args = list(),
    cluster = NULL, debug = FALSE) {

  # check the data are there.
  check.data(data)
  # check the number of bootstrap replicates.
  R = check.replicates(R)
  # check the size of each bootstrap sample.
  m = check.bootsize(m, data)
  # check the sim parameter.
  if (!(sim %in% c("ordinary", "parametric")))
    stop("the bootstrap simulation can be either 'ordinary' or 'parametric'.")
  # check debug.
  check.logical(debug)
  # check the learning algorithm.
  check.learning.algorithm(algorithm)
  # check the extra arguments for the learning algorithm.
  algorithm.args = check.learning.algorithm.args(algorithm.args)
  # check the custom statistic function.
  statistic = match.fun(statistic)
  # check the extra arguments for the statistic function.
  if (!is.list(statistic.args))
    statistic.args = as.list(statistic.args)

  # check the cluster.
  if (!is.null(cluster)) {

    check.cluster(cluster)

    # disable debugging, the slaves do not cat() here.
    if (debug) {

      warning("disabling debugging output in cluster-aware mode.")
      debug = FALSE

    }#THEN

  }#THEN

  bootstrap.backend(data = data, statistic = statistic, R = R, m = m,
    sim = sim, algorithm = algorithm, algorithm.args = algorithm.args,
    statistic.args = statistic.args, cluster = cluster, debug = debug)

}#BNBOOT

# compute arcs' strength via nonparametric bootstrap.
boot.strength = function(data, R = 200, m = nrow(data), algorithm,
    algorithm.args = list(), cpdag = TRUE, debug = FALSE) {

  # check the data are there.
  check.data(data)
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

  res = arc.strength.boot(data = data, R = R, m = m, algorithm = algorithm,
          algorithm.args = algorithm.args, arcs = NULL, cpdag = cpdag,
          debug = debug)

  # add extra information for strength.plot().
  res = structure(res, mode = "bootstrap", threshold = threshold(res),
          class = c("bn.strength", class(res)))

  return(res)

}#BOOT.STRENGTH

# perform cross-validation.
bn.cv = function(data, bn, loss = NULL, k = 10, algorithm.args = list(),
    loss.args = list(), fit = "mle", fit.args = list(), cluster = NULL,
    debug = FALSE) {

  # check the data are there.
  check.data(data)

  n = nrow(data)
  nodes = names(data)

  # check the number of splits.
  if (!is.positive.integer(k))
    stop("the number of splits must be a positive integer number.")
  if (k == 1)
    stop("the number of splits must be at least 2.")
  if (n < k)
    stop("insufficient sample size for ", k, " subsets.")

  # check the loss function.
  loss = check.loss(loss, data)
  # check the extra arguments passed down to the loss function.
  loss.args = check.loss.args(loss, bn, nodes, data, loss.args)
  # check the fitting method (maximum likelihood, bayesian, etc.)
  check.fitting.method(fit, data)
  # check the cluster.
  if (!is.null(cluster)) {

    check.cluster(cluster)

    # disable debugging, the slaves do not cat() here.
    if (debug) {

      warning("disabling debugging output in cluster-aware mode.")
      debug = FALSE

    }#THEN

  }#THEN

  if (is.character(bn)) {

    # check the learning algorithm.
    check.learning.algorithm(bn)
    # check whether it does return a DAG or not.
    if (!(bn %in% always.dag.result) && (loss == "pred"))
      stop(paste("this learning algorithm may result in a partially directed",
        "or undirected network, which is not handled by parameter fitting."))
    # check the extra arguments for the learning algorithm.
    algorithm.args = check.learning.algorithm.args(algorithm.args)

  }#THEN
  else if (is(bn, "bn")) {

    if (!identical(algorithm.args, list()))
      warning("no learning algorithm is used, so 'algorithm.args' will be ignored.")
    # check whether the data agree with the bayesian network.
    check.bn.vs.data(bn, data)
    # no parameters if the network structure is only partially directed.
    if (is.pdag(bn$arcs, nodes) && (loss == "pred"))
      stop("the graph is only partially directed.")
    # check bn.naive objects if any.
    if (is(bn, "bn.naive"))
      check.bn.naive(bn)

  }#THEN
  else {

    stop("bn must be the either label of a learning algorithm or a bn object.")

  }#ELSE

  crossvalidation(data = data, bn = bn, loss = loss, k = k,
    algorithm.args = algorithm.args, loss.args = loss.args,
    fit = fit, fit.args = fit.args, cluster = cluster, debug = debug)


}#BN.CV
