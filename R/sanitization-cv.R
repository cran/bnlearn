
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

  # check the number of splits.
  if (has.argument(method, "k", cv.extra.args))
    extra.args$k = check.cv.splits(extra.args$k, n = n)

  # check the number of runs.
  if (has.argument(method, "runs", cv.extra.args))
    extra.args$runs = check.cv.runs(extra.args$runs)

  # check the size of the test subsets in hold-out cross-validation.
  if (has.argument(method, "m", cv.extra.args))
    extra.args$m = check.cv.holdout.size(extra.args$m, n = n)

  # warn about unused arguments (before adding more below).
  check.unused.args(extra.args, cv.extra.args[[method]])

  # check custom folds.
  if (has.argument(method, "folds", cv.extra.args)) {

    extra.args$folds = check.cv.folds(extra.args$folds, n = n)
    extra.args$runs = length(extra.args$folds)

  }#THEN

  return(extra.args)

}#CHECK.CV.ARGS

# check the number of splits.
check.cv.splits = function(k, n) {

  # the default is 10 folds, or leave-one-out cross-validation if there are
  # fewer that 10 observations.
  if (is.null(k))
    k = min(10, n)

  if (!is.positive.integer(k))
    stop("the number of splits must be a positive integer number.")
  if (k == 1)
    stop("the number of splits must be at least 2.")
  if (n < k)
    stop("insufficient sample size for ", k, " subsets.")

  return(k)

}#CHECK.CV.SPLITS

# check the number of runs.
check.cv.runs = function(runs)  {

  if (is.null(runs))
    return(1)

  if (!is.positive.integer(runs))
    stop("the number of runs must be a positive integer number.")

  return(runs)

}#CHECK.CV.RUNS

# check the size of the test subsets in hold-out cross-validation.
check.cv.holdout.size = function(m, n) {

  if (is.null(m))
    return(ceiling(n / 10))

  if (!is.positive.integer(m))
    stop("the size of the test subset must be a positive integer number.")
  if (m >= n)
    stop("insufficient sample size for a test subset of size ", m, ".")

  return(m)

}#CHECK.CV.HOLDOUT.SIZE

# check custom folds.
check.cv.folds = function(folds, n) {

  if (is.null(folds))
    stop("custom folds are missing, with no default.")

  check.single.run = function(folds, n) {

    if (!is.list(folds))
      stop("folds must be specified via a list of indices.")
    if (length(folds) < 2)
      stop("at least two folds are needed.")
    if (any(sapply(folds, length) == 0))
      stop("some folds contain no observations.")
    if (any(!sapply(folds, is.positive.vector)))
      stop("observation indices must be positive integer numbers.")

    merged = unlist(folds)

    if (any(duplicated(merged)))
      stop("some observations are included in more than one fold.")
    if (any(merged > n))
      stop("observation indices are too high (sample size is ", n, ").")
    if (length(merged) != n)
      stop("not all observations are assigned to a fold.")

    return(folds)

  }#CHECK.SINGLE.RUN

  if (!is.list(folds))
    stop("folds must be specified via a list of indices.")

  if (all(sapply(folds, is.list))) {

    lapply(folds, check.single.run, n = n)
    return(folds)

  }#THEN
  else {

    check.single.run(folds = folds, n = n)
    return(list(folds))

  }#ELSE

}#CHECK.CV.FOLDS

