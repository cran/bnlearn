
# descriptive statistics of network's variability.
bn.var = function(x, method) {

  # match the test statistic.
  method.string = method
  method = check.mvber.vartest(method)

  if (is(x, "mvber.moments")) {

    mvber.var.backend(x$covariance, method = method)

  }#THEN
  else {

    # check the covariance matrix.
    check.covariance(x)

    mvber.var.backend(x, method = method)

  }#ELSE

}#BN.VAR

# estimate the first two moments of the  multivariate Bernoulli distribution
# associated with a set of bootstrapped Bayesian networks.
bn.moments = function(data, R = 200, m = nrow(data), algorithm,
    algorithm.args = list(), reduce = NULL, debug = FALSE) {

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

  res = mvber.moments.backend(data = data, R = R, m = m, algorithm = algorithm,
          algorithm.args = algorithm.args, arcs = NULL, debug)

  # reduce the return value either removing any arc which did not appear
  # in any bootstrap sample ("first") or that has zero variance ("second").
  if (!is.null(reduce)) {

    if (!is.string(reduce) || !(reduce %in% c("first", "second")))
      stop("unknown criterion for matrix reduction.")

    if (reduce == "first")
      kill = which(res$expected == 0)
    else if (reduce == "second")
      kill = which(diag(res$covariance) == 0)

    if (length(kill) > 0) {

      res$expected = res$expected[-kill]
      res$covariance = res$covariance[-kill, -kill, drop = FALSE]

    }#THEN

  }#THEN

  return(structure(res, class = c("mvber.moments", class(res)), R = R, m = m))

}#BN.MOMENTS

