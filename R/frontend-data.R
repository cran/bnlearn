
# transform continuous data into discrete ones.
discretize = function(data, method, breaks = 3, ordered = FALSE, ...,
    debug = FALSE) {

  # check the label of the discretization method.
  method = check.discretization.method(method)
  # general check on the data.
  data = check.data(data, allow.missing = TRUE, stop.if.all.missing = TRUE)

  # the data should include at least some continuous variables, otherwise we
  # have nothing to do.
  if (attr(data, "metadata")$type %in% discrete.data.types) {

    # ensure that the attribute with the metadata set by check.data() is removed.
    attr(data, "metadata") = NULL

    warning("at least one variable should be continuous")
    return(data)

  }#THEN

  # check the number of breaks.
  if (length(breaks) == 1)
    breaks = rep(breaks, ncol(data))
  else if (length(breaks) != ncol(data))
    stop("the 'breaks' vector must have an element for each variable in the data.")
  if (!is.positive.vector(breaks))
    stop("the number(s) of breaks must be positive integer number(s).")
  if (any(breaks == 1))
    stop("the return value must have at least two levels for each variable.")

  # check whether/which discretized variables should be ordered factors.
  if (length(ordered) == 1)
    ordered = rep(ordered, ncol(data))
  else if (length(ordered) != ncol(data))
    stop("the 'ordered' vector must have an element for each variable in the data.")
  if (!is.logical.vector(ordered))
    stop("the elements of the 'ordered' vector be logical values.")

  # check the data.
  if (method == "hartemink") {

    # check that the data contains at least two columns, otherwise there is
    # nothing to compute mutual information from.
    if (ncol(data) < 2)
      stop("at least two variables are needed to compute mutual information.")

  }#THEN

  # check debug.
  check.logical(debug)
  # check the extra arguments.
  extra.args = check.discretization.args(method, data, breaks, list(...))

  discretize.backend(data = data, method = method, breaks = breaks,
    ordered = ordered, extra.args = extra.args, debug = debug)

}#DISCRETIZE

# screen the data for highly correlated variables.
dedup = function(data, threshold = 0.90, debug = FALSE) {

  # check the data (only continuous data are supported).
  data = check.data(data, allowed.types = continuous.data.types,
           allow.missing = TRUE)
  # check the correlation threshold.
  if (missing(threshold))
    threshold = 0.90
  else if (!is.probability(threshold))
    stop("the correlation threshold must be a number between 0 and 1.")
  # check debug.
  check.logical(debug)

  dedup.backend(data = data, threshold = threshold, debug = debug)

}#DEDUP

# configurations of sets of discrete variables.
configs = function(data, all = TRUE) {

  # check the data (only discrete data are supported).
  data = check.data(data, allowed.types = discrete.data.types,
           allow.missing = TRUE, allow.levels = TRUE)
  # check the "all configurations" flag.
  check.logical(all)

  configurations(data, factor = TRUE, all = all)

}#CONFIGS
