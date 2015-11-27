
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
  if (sim %!in% c("ordinary", "parametric"))
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

      warning("disabling debugging output for parallel computing.")
      debug = FALSE

    }#THEN

  }#THEN

  bootstrap.backend(data = data, statistic = statistic, R = R, m = m,
    sim = sim, algorithm = algorithm, algorithm.args = algorithm.args,
    statistic.args = statistic.args, cluster = cluster, debug = debug)

}#BNBOOT

# compute arcs' strength via nonparametric bootstrap.
boot.strength = function(data, cluster = NULL, R = 200, m = nrow(data),
    algorithm, algorithm.args = list(), cpdag = TRUE, debug = FALSE) {

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

  # check the cluster.
  if (!is.null(cluster)) {

    check.cluster(cluster)

    # enter in cluster-aware mode.
    cluster.aware = TRUE
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
  res = structure(res, mode = "bootstrap", threshold = threshold(res),
          class = c("bn.strength", class(res)))

  return(res)

}#BOOT.STRENGTH

# perform cross-validation.
bn.cv = function(data, bn, loss = NULL, k = 10, m, runs = 1, 
    algorithm.args = list(), loss.args = list(), fit = "mle",
    fit.args = list(), method = "k-fold", cluster = NULL, debug = FALSE) {

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
  # check the number of runs.
  if (!is.positive.integer(runs))
    stop("the number of runs must be a positive integer number.")

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
  # check the cross-validation method.
  if (!is.string(method) && (method %!in% c("k-fold", "hold-out") ))
    stop("valid cross-validation methods are 'k-fold' and 'hold-out'.")
  # check the size of the test subsets in hold-put cross-validation.
  if (method == "hold-out") {

    if (missing(m))
      m = ceiling(n / 10)
    if (!is.positive.integer(m))
      stop("the size of the test subset must be a positive integer number.")
    if (m >= n)
      stop("insufficient sample size for a test subset of size ", m, ".")

  }#THEN
  else {

    if (!missing(m))
      warning("'m' is only required for hold-out cross-validation, ignoring.")

  }#ELSE

  if (is.character(bn)) {

    # check the learning algorithm.
    check.learning.algorithm(bn)
    # check the loss function.
    loss = check.loss(loss, data, bn)
    # check whether it does return a DAG or not.
    if ((bn %!in% always.dag.result) && (loss == "pred"))
      stop("this learning algorithm may result in a partially directed",
        " or undirected network, which is not handled by parameter fitting.")
    # check the extra arguments for the learning algorithm.
    algorithm.args = check.learning.algorithm.args(algorithm.args)
    # since we have no other way to guess, copy the label of the target
    # variable from the parameters of the classifier.
    if ((loss == "pred") && (is.null(loss.args$target)) && (bn %in% classifiers))
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
    # no parameters if the network structure is only partially directed.
    if (is.pdag(bn$arcs, nodes) && (loss == "pred"))
      stop("the graph is only partially directed.")
    # check bn.naive objects if any.
    if (is(bn, "bn.naive"))
      check.bn.naive(bn)
    # check the extra arguments passed down to the loss function.
    loss.args = check.loss.args(loss, bn, nodes, data, loss.args)

  }#THEN
  else {

    stop("bn must be either the label of a learning algorithm or a bn object.")

  }#ELSE

  # allocate and populate the return value.
  result = structure(vector(runs, mode = "list"), class = "bn.kcv.list")

  for (r in seq(runs))
    result[[r]] = crossvalidation(data = data, bn = bn, loss = loss, k = k,
                    m = m, algorithm.args = algorithm.args,
                    loss.args = loss.args, fit = fit, fit.args = fit.args,
                    method = method, cluster = cluster, debug = debug)

  # return a bn.kcv object (for a single run) or a bn.kcv.list object (for
  # multiple runs).
  if (runs == 1)
    return(result[[1]])
  else
    return(result)

}#BN.CV
