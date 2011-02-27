
# transform continuous data into discrete ones.
discretize = function(x, method, breaks = 3, ..., debug = FALSE) {

  # check the label of the discretization method.
  method = check.discretization.method(method)
  # check the data.
  if (method %in% c("quantile", "interval")) {

    check.data(x, allow.mixed = TRUE)
    # check the number of breaks.
    if (length(breaks) == 1) {

      # get an array of the correct length.
      breaks = rep(breaks, ncol(data))

    }#THEN
    else if (length(breaks) != ncol(data))
      stop("the 'breaks' vector must have an element for each variable in the data.")
    if (!is.positive.vector(breaks))
      stop("the numbers of breaks must be positive integer numbers.")

  }#THEN
  else if (method == "hartemink") {

    check.data(x, allow.mixed = FALSE)
    # Haretmink's method requires all the variables to be continuous.
    if (!is.data.continuous(x))
      stop("all variables must be continuous.")
    # check the number of breaks.
    if (!is.positive.integer(breaks))
      stop("the number of breaks must be a positive integer number.")

  }#THEN

  # check debug.
  check.logical(debug)
  # check the extra arguments.
  extra.args = check.discretization.args(method, x, list(...))

  discretize.backend(data = x, method = method, breaks = breaks, 
    extra.args = extra.args, debug = debug)

}#DISCRETIZE

