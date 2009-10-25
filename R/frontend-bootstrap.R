
# generic frontend to {non,}parametric bootstrap.
bnboot = function(data, statistic, R = 200, m = nrow(data),
    sim = "ordinary", algorithm, algorithm.args = list(),
    statistic.args = list(), debug = FALSE) {

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

  bootstrap.backend(data = data, statistic = statistic, R = R, m = m,
    sim = sim, algorithm = algorithm, algorithm.args = algorithm.args,
    statistic.args = statistic.args, debug = debug)

}#BNBOOT

boot.strength = function(data, R = 200, m = nrow(data),
    algorithm, algorithm.args = list(), debug = FALSE) {

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

  arc.strength.boot(data = data, R = R, m = m, algorithm = algorithm,
    algorithm.args = algorithm.args, arcs = NULL, debug = debug)

}#BOOT.STRENGTH

