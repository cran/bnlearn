# check bootstrap arguments (when they are passed as variable length args).
check.bootstrap.args = function(extra.args, network, data) {

  # check the number of bootstrap replicates.
  extra.args$R = check.replicates(extra.args$R)
  # check the size of each bootstrap sample.
  extra.args$m = check.bootsize(extra.args$m, data)
  # check the learning algorithm.
  algorithm = check.learning.algorithm(extra.args[["algorithm"]], bn = network)
  # check the extra arguments for the learning algorithm.
  algorithm.args = check.learning.algorithm.args(extra.args[["algorithm.args"]],
                     algorithm = algorithm, bn = network)

  extra.args[["algorithm"]] = algorithm
  extra.args[["algorithm.args"]] = algorithm.args

  # remap additional arguments used in hybrid algorithms.
  if (algorithm %in% hybrid.algorithms) {

    # there's no need to sanitize these parameters, it's done either in
    # bnlearn() or in greedy.search() already.
    if (is.null(extra.args[["algorithm.args"]]$restrict))
      extra.args[["algorithm.args"]]$restrict = network$learning$restrict
    if (is.null(extra.args[["algorithm.args"]]$maximize))
      extra.args[["algorithm.args"]]$maximize = network$learning$maximize
    if (is.null(extra.args[["algorithm.args"]]$test))
      extra.args[["algorithm.args"]]$test = network$learning$rstest
    if (is.null(extra.args[["algorithm.args"]]$score))
      extra.args[["algorithm.args"]]$score = network$learning$maxscore

  }#THEN

  # warn about unused arguments.
  check.unused.args(extra.args, c("R", "m", "algorithm", "algorithm.args"))

  return(extra.args)

}#CHECK.BOOTSTRAP.ARGS

# check the number of bootstrap replicates.
check.replicates = function(R, default = 200) {

  if (missing(R) || is.null(R))
    R = default
  else if (!is.positive.integer(R))
    stop("the number of bootstrap replicates must be a positive integer.")

  return(R)

}#CHECK.REiPLICATES

# check the size of bootstrap replicates.
check.bootsize = function(m, data) {

  if (missing(m) || is.null(m))
    m = nrow(data)
  else if (!is.positive.integer(m))
    stop("bootstrap sample size must be a positive integer.")

  return(m)

}#CHECK.BOOTSIZE
