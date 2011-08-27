
# transform continuous data into discrete ones.
discretize = function(x, method, breaks = 3, ..., debug = FALSE) {

  # check the label of the discretization method.
  method = check.discretization.method(method)
  # check the data.
  if (method %in% c("quantile", "interval")) {

    # general check on the data.
    check.data(x, allow.mixed = TRUE)
    # check the number of breaks.
    if (length(breaks) == 1) {

      # get an array of the correct length.
      breaks = rep(breaks, ncol(x))

    }#THEN
    else if (length(breaks) != ncol(x))
      stop("the 'breaks' vector must have an element for each variable in the data.")
    if (!is.positive.vector(breaks))
      stop("the numbers of breaks must be positive integer numbers.")

  }#THEN
  else if (method == "hartemink") {

    # general check on the data.
    check.data(x, allow.mixed = FALSE)

    # check the number of breaks.
    if (!is.positive.integer(breaks))
      stop("the number of breaks must be a positive integer number.")
    if (breaks == 1)
      stop("the return value must have at least two levels for each variable.")

    # check the data types.
    if (is.data.discrete(x)) {

      nlvls = unique(sapply(x, nlevels))
      # this is implicit in Hartemink's definition of the algorithm.
      if (length(nlvls) > 1)
        stop("all variables must have the same number of levels.")
      # we only aggregate levels, so must have enough of them.
      if (nlvls < breaks)
        stop("too many breaks, at most ", nlvls, " required.")

    }#THEN
    else if (!is.data.continuous(x)) {

      stop("all variables must be continuous.")

    }#THEN

  }#THEN

  # check debug.
  check.logical(debug)
  # check the extra arguments.
  extra.args = check.discretization.args(method, x, breaks, list(...))

  discretize.backend(data = x, method = method, breaks = breaks,
    extra.args = extra.args, debug = debug)

}#DISCRETIZE

