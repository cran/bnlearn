
# check the label of a learning algorithm.
check.learning.algorithm = function(algorithm, class = "all", bn) {

  # select the right class of algorithms.
  ok = character(0)

  if ("constraint" %in% class)
    ok = c(ok, constraint.based.algorithms)
  if ("markov.blanket" %in% class)
    ok = c(ok, markov.blanket.algorithms)
  if ("neighbours" %in% class)
    ok = c(ok, local.search.algorithms)
  if ("score" %in% class)
    ok = c(ok, score.based.algorithms)
  if ("mim" %in% class)
    ok = c(ok, mim.based.algorithms)
  if ("classifier" %in% class)
    ok = c(ok, classifiers)
  if ("all" %in% class)
    ok = available.learning.algorithms

  if (missing(algorithm) || is.null(algorithm)) {

    # use the one specified by the bn object as the default.
    if (missing(bn))
      stop("the learning algorithm must be a character string.")
    else if (is(bn, "bn"))
      algorithm = bn$learning$algo

  }#THEN
  else {

    check.label(algorithm, choices = ok, labels = method.labels,
      argname = "learning algorithm", see = "bnlearn-package")

  }#ELSE

  return(algorithm)

}#CHECK.LEARNING.ALGORITHM

# check size of the largest conditioning set in the independence tests.
check.largest.sx.set = function(max.sx, data) {

  if (missing(max.sx) || is.null(max.sx)) {

    max.sx = Inf

  }#THEN
  else if (!isTRUE(all.equal(max.sx, Inf))) {

    if (!is.positive.integer(max.sx))
      stop("max.sx must be a positive integer number.")
    if (max.sx >= ncol(data))
      warning("max.sx should be lower than the number of nodes, the limit will be ignored.")

  }#ELSE

  return(max.sx)

}#CHECK.LARGEST.SX.SET

# check the threshold on the maximum number of parents.
check.maxp = function(maxp, data) {

  if (missing(maxp) || is.null(maxp)) {

    maxp = Inf

  }#THEN
  else if (!isTRUE(all.equal(maxp, Inf))) {

    if (!is.positive.integer(maxp))
      stop("maxp must be a positive integer number.")
    if (maxp >= ncol(data))
      warning("maximum number of parents should be lower than the number of nodes, the limit will be ignored.")

  }#ELSE

  return(as.numeric(maxp))

}#CHECK.MAXP

# check parameters related to the random restart functions.
check.restart = function(restart) {

  # set the default value if not specified.
  if (missing(restart) || is.null(restart) || (restart == 0))
      return(0)

  if (!is.positive.integer(restart))
    stop("the number of random restarts must be a non-negative numeric value.")
  else
    return(restart)

}#CHECK.RESTART

check.perturb = function(perturb) {

  # set the default value if not specified.
  if (missing(perturb) || is.null(perturb))
      return(1)

  if (!is.positive.integer(perturb))
    stop("the number of changes at each radom restart must be a non-negative numeric value.")
  else
    return(perturb)

}#CHECK.PERTURB

# check the maximum number of iterations.
check.max.iter = function(max.iter) {

  # set the default value if not specified.
  if (missing(max.iter) || is.null(max.iter))
    return(Inf)

  if ((max.iter != Inf) && !is.positive.integer(max.iter))
    stop("the maximum number of iterations must be a positive integer number.")
  else
    return(max.iter)

}#CHECK.MAX.ITER

# check arguments related to the tabu list.
check.tabu = function(tabu) {

  # set the default value if not specified.
  if (missing(tabu) || is.null(tabu))
    return(10)

  if (!is.positive.integer(tabu))
    stop("the length of the tabu list must be a positive integer number.")
  else
    return(tabu)

}#CHECK.TABU

check.max.tabu = function(max, tabu) {

  if (missing(max) || is.null(max))
    return(tabu)

  # check the number of iterations the algorithm can perform without
  # improving the best network score.
  if (!is.positive.integer(max))
    stop("the maximum number of iterations without any score improvement must be a positive integer number.")
  # the tabu list should be longer than that, otherwise the search can do a
  # U-turn and return to the local maximum it left before (thus creating an
  # endless loop).
  if (max > tabu)
    stop("the maximum number of iterations without any score improvement must not be grater than the length of the tabu list.")

  return(max)

}#CHECK.MAX.TABU

# check the arguments of a learning algorithm (for use in bootstrap).
check.learning.algorithm.args = function(args, algorithm, bn) {

  # convert args into a list, if it's not one already.
  if (!is.list(args))
      args = as.list(args)

  # if a reference bn is specified, guess as many parameters as possbile.
  if (!(missing(algorithm) || missing(bn))) {

    # use the same score/conditional independence test.
    if (algorithm %in% constraint.based.algorithms) {

      # it's essential to check it's actually an independence test,
      # it could be a score function or NA.
      if ("test" %!in% names(args))
        if (bn$learning$test %in% available.tests)
          args$test = bn$learning$test

      # set the appropriate value for the optimization flag.
      if ("optimized" %!in% names(args))
        args$optimized = bn$learning$optimized

      # pass along all the parameters in bn$learning$args.
      if (length(bn$learning$args) > 0) {

        if ("alpha" %!in% names(args))
          args$alpha = bn$learning$args$alpha

        if ("test" %in% names(args) && ("B" %!in% names(args)))
          if (args$test %in% resampling.tests)
            args$B = bn$learning$args$B

      }#THEN

    }#THEN
    else if (algorithm %in% score.based.algorithms) {

      if ("score" %!in% names(args))
        if (bn$learning$test %in% available.scores)
          args$score = bn$learning$test

      # set the appropriate value for the optimization flag.
      if ("optimized" %!in% names(args))
        args$optimized = bn$learning$optimized

      # pass along the relevant parameters in bn$learning$args if the score
      # function is the same (hint: different scores have paramenters with
      # the same name but different meanings).
      if (("score" %in% names(args)) && (args$score == bn$learning$test))
        for (arg in names(bn$learning$args))
          if ((arg %!in% names(args)) && (arg %in% (score.extra.args[[args$score]])))
            args[[arg]] = bn$learning$args[[arg]]

    }#THEN

    # pass along whitelist and blacklist.
    if (!is.null(bn$learning$whitelist))
      args$whitelist = bn$learning$whitelist
    if (!is.null(bn$learning$blacklist))
      args$blacklist = bn$learning$blacklist

  }#THEN

  # remove any spurious x arguments, the data are provided by the bootstrap.
  if ("x" %in% names(args)) {

    args$x = NULL

    warning("removing 'x' from 'algorithm.args', the data set is provided by the bootstrap sampling.")

  }#THEN

  return(args)

}#CHECK.LEARNING.ALGORITHM.ARGS

