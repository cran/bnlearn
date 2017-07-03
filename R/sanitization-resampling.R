
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
check.bootsize = function(m, data, default = nrow(data)) {

  if (missing(m) || is.null(m))
    m = default
  else if (!is.positive.integer(m))
    stop("bootstrap sample size must be a positive integer.")

  return(m)

}#CHECK.BOOTSIZE

# check the method used for cross-validation.
check.cv.method = function(method) {

  if (!missing(method) && !is.null(method)) {

    check.label(method, choices = available.cv.methods,
      labels = cv.labels, argname = "cross-validation method", see = "bn.cv")

    return(method)

  }#THEN
  else {

    return("k-fold")

  }#ELSE

}#CHECK.CV.METHOD

# sanitize the extra arguments passed to cross-validation methods.
check.cv.args = function(method, extra.args, data) {

  n = nrow(data)

  if (method %in% c("k-fold", "hold-out")) {

    # check the number of splits.
    if (!is.null(extra.args$k)) {

      if (!is.positive.integer(extra.args$k))
        stop("the number of splits must be a positive integer number.")
      if (extra.args$k == 1)
        stop("the number of splits must be at least 2.")
      if (n < extra.args$k)
        stop("insufficient sample size for ", extra.args$k, " subsets.")

    }#THEN
    else {

      extra.args$k = 10

    }#ELSE

    # check the number of runs.
    if (!is.null(extra.args$runs)) {

      if (!is.positive.integer(extra.args$runs))
        stop("the number of runs must be a positive integer number.")

    }#TTHEN
    else {

      extra.args$runs = 1

    }#ELSE

  }#THEN

  if (method == "hold-out") {

    # check the size of the test subsets in hold-put cross-validation.
    if (!is.null(extra.args$m)) {

      if (!is.positive.integer(extra.args$m))
        stop("the size of the test subset must be a positive integer number.")
      if (extra.args$m >= n)
        stop("insufficient sample size for a test subset of size ",
          extra.args$m, ".")

    }#THEN
    else {

      extra.args$m = ceiling(n / 10)

    }#ELSE

  }#THEN

  if (method == "custom-folds") {

    if (!is.null(extra.args$folds)) {

      if (!is.list(extra.args$folds))
        stop("folds must be specified via a list of indices.")
      if (length(extra.args$folds) < 2)
        stop("at least two folds are needed.")
      if (any(sapply(extra.args$folds, length) == 0))
        stop("some folds contain no observations.")
      if (any(!sapply(extra.args$folds, is.positive.vector)))
        stop("observation indices must be positive integer numbers.")

      merged = unlist(extra.args$folds)

      if (any(duplicated(merged)))
        stop("some observations are included in more than one fold.")
      if (any(merged > n))
        stop("observation indices are too high (sample size is ", n, ").")
      if (length(merged) != n)
        stop("not all observations are assigned to a fold.")

    }#THEN
    else {

      stop("custom folds are missing, with no default.")

    }#ELSE

  }#THEN

  # warn about unused arguments.
  check.unused.args(extra.args, cv.extra.args[[method]])

  return(extra.args)

}#CHECK.CV.ARGS

