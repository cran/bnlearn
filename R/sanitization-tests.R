
# check test labels.
check.test = function(test, data) {

  # check which type of data we are dealing with.
  type = data.type(data)

  if (!missing(test) && !is.null(test)) {

    # check the test label.
    check.label(test, choices = available.tests, labels = test.labels,
      argname = "conditional independence test", see = "bnlearn-package")
    # check if it's the right test for the data (discrete, continuous).
    if ((type != "ordered") && (test %in% available.ordinal.tests))
      stop("test '", test, "' may be used with ordinal data only.")
    if ((type %!in% discrete.data.types) && (test %in% available.discrete.tests))
      stop("test '", test, "' may be used with discrete data only.")
    if ((type != "continuous") && (test %in% available.continuous.tests))
      stop("test '", test, "' may be used with continuous data only.")
    if ((type != "mixed-cg") && (test %in% available.mixedcg.tests))
      stop("test '", test, "' may be used with a mixture of continuous and discrete data only.")

    return(test)

  }#THEN
  else {

    if (type == "ordered")
      return("jt")
    else if (type %in% c("factor", "mixed-do"))
      return("mi")
    else if (type == "continuous")
      return("cor")
    else if (type == "mixed-cg")
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

