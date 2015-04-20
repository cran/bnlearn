
# transform continuous data into discrete ones.
discretize = function(data, method, breaks = 3, ordered = FALSE, ..., debug = FALSE) {

  # check the label of the discretization method.
  method = check.discretization.method(method)
  # general check on the data.
  type = check.data(data)

  # check the data.
  if (method %in% c("quantile", "interval")) {

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

    # check the number of breaks.
    if (!is.positive.integer(breaks))
      stop("the number of breaks must be a positive integer number.")
    if (breaks == 1)
      stop("the return value must have at least two levels for each variable.")

    # check the data types.
    if (type %in% discrete.data.types) {

      nlvls = unique(sapply(data, nlevels))
      # this is implicit in Hartemink's definition of the algorithm.
      if (length(nlvls) > 1)
        stop("all variables must have the same number of levels.")
      # we only aggregate levels, so must have enough of them.
      if (nlvls < breaks)
        stop("too many breaks, at most ", nlvls, " required.")

    }#THEN
    else if (type != "continuous") {

      stop("all variables must be continuous.")

    }#THEN

  }#THEN

  # check debug and ordered.
  check.logical(debug)
  check.logical(ordered)
  # check the extra arguments.
  extra.args = check.discretization.args(method, data, breaks, list(...))

  discretize.backend(data = data, method = method, breaks = breaks,
    ordered = ordered, extra.args = extra.args, debug = debug)

}#DISCRETIZE

# Pena's relevant nodes feature selection.
relevant = function(target, context, data, test, alpha, B, debug = FALSE) {

  # check the data.
  check.data(data, allowed.types = c(discrete.data.types, continuous.data.types))
  # a valid node is needed.
  check.nodes(nodes = target, graph = data, max.nodes = 1)
  # an optional valid node is needed.
  if (!missing(context)) {

    check.nodes(nodes = context, graph = data)

    if (length(intersect(target, context)) > 0)
      stop("target and context nodes must be disjoint sets.")

  }#THEN
  else {

    context = NULL

  }#ELSE
  # check the test label.
  test = check.test(test, data)
  # check B (the number of eprmutations).
  B = check.B(B, test)
  # check alpha.
  alpha = check.alpha(alpha)
  # check debug.
  check.logical(debug)

  pena.backend(target = target, context = context, data = data, test = test,
    alpha = alpha, B = B, debug = debug)

}#RELEVANT

# screen the data for highly correlated variables.
dedup = function(data, threshold = 0.90, debug = FALSE) {

  # check the data (only continuous data are supported).
  type = check.data(data, allowed.types = continuous.data.types)
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
  check.data(data, allowed.types = discrete.data.types)
  # check the "all configurations" flag.
  check.logical(all)

  configurations(data, factor = TRUE, all = all)

}#CONFIGS
