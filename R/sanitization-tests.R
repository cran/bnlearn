
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

# check the the target nominal type I error rate
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

# check the number of permutation/boostrap samples.
check.B = function(B, criterion) {

  if (criterion %in% resampling.tests) {

    if (!missing(B) && !is.null(B)) {

      if (!is.positive.integer(B))
        stop("the number of permutations/bootstrap replications must be a positive integer number.")

      B = as.integer(B)

    }#THEN
    else {

      if (criterion %in% semiparametric.tests)
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

