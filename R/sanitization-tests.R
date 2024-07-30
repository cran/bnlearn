
# check test labels.
check.test = function(test, data) {

  data.type = attr(data, "metadata")$type
  if (is.null(data.type))
    data.type = data.type(data)

  if (!missing(test) && !is.null(test)) {

    # check the test label.
    check.label(test, choices = available.tests, labels = test.labels,
      argname = "conditional independence test", see = "bnlearn-package")
    # check if it's the right test for the data (discrete, continuous).
    if ((data.type != "ordered") &&
        (test %in% available.ordinal.tests))
      stop("test '", test, "' may only be used with ordinal data.")
    if ((data.type %!in% discrete.data.types) &&
        (test %in% available.discrete.tests))
      stop("test '", test, "' may only be used with discrete data.")
    if ((data.type != "continuous") &&
        (test %in% available.continuous.tests))
      stop("test '", test, "' may only be used with continuous data.")
    if ((data.type != "mixed-cg") &&
        (test %in% available.mixedcg.tests))
      stop("test '", test, "' may only be used with a mixture of continuous and discrete data.")

    return(test)

  }#THEN
  else {

    if (data.type == "ordered")
      return("jt")
    else if (data.type %in% c("factor", "mixed-do"))
      return("mi")
    else if (data.type == "continuous")
      return("cor")
    else if (data.type == "mixed-cg")
      return("mi-cg")

  }#ELSE

}#CHECK.TEST

# check the the target nominal type I error rate.
check.alpha = function(alpha, network = NULL) {

  # check the the target nominal type I error rate
  if (!missing(alpha) && !is.null(alpha)) {

    # validate alpha.
    if (!is.probability(alpha))
      stop("alpha must be a numerical value in [0,1].")

  }#THEN
  else {

    # check if there is an alpha value stored in the bn object;
    # otherwise use the usual 0.05 value.
    if (!is.null(network$learning$args$alpha))
      alpha = network$learning$args$alpha
    else
      alpha = 0.05

  }#ELSE

  return(alpha)

}#CHECK.ALPHA

# sanitize the extra arguments passed to the conditional independence tests.
check.test.args = function(test, network = NULL, data, extra.args) {

  # check the imaginary sample size.
  if (has.argument(test, "B", test.extra.args))
    extra.args[["B"]] =
      check.B(B = extra.args[["B"]], network = network, criterion = test)

  # check the R function implementing a custom test.
  if (has.argument(test, "fun", test.extra.args))
    extra.args[["fun"]] =
      check.custom.test.function(fun = extra.args[["fun"]], network = network)
  if (has.argument(test, "args", test.extra.args))
    extra.args[["args"]] =
      check.custom.test.arguments(args = extra.args[["args"]], network = network)

  # warn about and remove unused arguments.
  extra.args = check.unused.args(extra.args, test.extra.args[[test]])

  return(extra.args)

}#CHECK.TEST.ARGS

# check the number of permutation/boostrap samples.
check.B = function(B, network, criterion) {

  if (criterion %in% resampling.tests) {

    if (!missing(B) && !is.null(B)) {

      if (!is.positive.integer(B))
        stop("the number of permutations/bootstrap replications must be a positive integer number.")

      B = as.integer(B)

    }#THEN
    else {

      if (!is.null(network$learning$args$B))
        B = network$learning$args$B
      else if (criterion %in% semiparametric.tests)
        B = 100L
      else
        B = 5000L

    }#ELSE

  }#THEN
  else {

    if (!missing(B) && !is.null(B))
      warning("this test does not require any permutations/bootstrap resampling, ignoring B.\n")

    B = NULL

  }#ELSE

  return(B)

}#CHECK.B

# check the R function implementing a custom score.
check.custom.test.function = function(fun, network) {

  # there is no possible default value.
  if (is.null(fun)) {

    if (!is.null(network$learning$args$fun))
      fun = network$learning$args$fun
    else
      stop("missing the custom test function.")

  }#THEN
  else {

    # check the argument list.
    fun.arguments = names(formals(fun))
    if (!identical(fun.arguments, c("x", "y", "z", "data", "args")))
      stop("the custom test function must have signature function(x, y, z, data, args).")

  }#ELSE

  return(fun)

}#CHECK.CUSTOM.TEST.FUNCTION

# check the additional argument list passed to a custom test.
check.custom.test.arguments = function(args, network) {

  # default to an empty argument list.
  if (is.null(args)) {

    if (!is.null(network$learning$args$args))
      args = network$learning$args$args
    else
      args = list()

  }#THEN
  else {

    if (!is.list(args))
      stop("the arguments for the custom test must be passed as a list.")

  }#ELSE

  return(args)

}#CHECK.CUSTOM.TEST.ARGUMENTS
