
# check the method used to discretize the data.
check.discretization.method = function(method) {

  if (!missing(method) && !is.null(method)) {

    check.label(method, choices = available.discretization.methods,
      labels = discretization.labels, argname = "discretization method",
      see = "discretize")

    return(method)

  }#THEN
  else {

    return("quantile")

  }#ELSE

}#CHECK.DISCRETIZATION.METHOD

# sanitize the extra arguments passed to discretization methods.
check.discretization.args = function(method, data, breaks, extra.args) {

  # check which type of data we are dealing with.
  type = data.type(data)

  # check the initial discretization algorithm.
  if (has.argument(method, "idisc", discretization.extra.args))
    extra.args[["idisc"]] = check.idisc(extra.args[["idisc"]])

  # check the initial number of breaks.
  if (has.argument(method, "ibreaks", discretization.extra.args))
    extra.args[["ibreaks"]] =
      check.ibreaks(extra.args[["ibreaks"]], breaks = breaks, data = data)

  # warn about and remove unused arguments.
  extra.args = check.unused.args(extra.args, discretization.extra.args[[method]])

  return(extra.args)

}#CHECK.DISCRETIZATION.ARGS

# check the initial discretization algorithm.
check.idisc = function(idisc) {

  if (!is.null(idisc)) {

    other.methods = setdiff(available.discretization.methods, "hartemink")

    check.label(idisc, choices = other.methods, labels = discretization.labels,
      argname = "initial discretization method", see = "discretize")

  }#THEN
  else {

    # default to quantile discretization as per Hartemink's recommendation.
    idisc = "quantile"

  }#ELSE

  return(idisc)

}#CHECK.IDISC

# check the initial number of breaks.
check.ibreaks = function(ibreaks, breaks, data) {

  if (!is.null(ibreaks)) {

    if (length(ibreaks) == 1)
      ibreaks = rep(ibreaks, ncol(data))
    else if (length(ibreaks) != ncol(data))
      stop("the 'ibreaks' vector must have an element for each variable in the data.")
    if (!is.positive.vector(ibreaks))
      stop("the number of initial breaks must be a (vector of) positive integer number(s).")
    if (any(ibreaks < breaks))
      stop("insufficient number of initial breaks for variables ",
           paste0(names(data)[ibreaks < breaks], collapse = ", "), ".")
    else if (any(ibreaks == breaks))
      warning("the initial number of breaks is identical to the final number of breaks.")

  }#THEN
  else {

    ndata = nrow(data)

    if (ndata > 500)
      ibreaks = 50
    else if (ndata > 100)
      ibreaks = 20
    else if (ndata > 50)
      ibreaks = 10
    else if (ndata > 10)
      ibreaks = 5
    else
      ibreaks = ndata

    ibreaks = rep(ibreaks, ncol(data))

  }#ELSE

  return(ibreaks)

}#CHECK.IBREAKS

# check the method used to deduplicate the data.
check.deduplication.method = function(method, data) {

  metadata = attr(data, "metadata")

  if (!missing(method) && !is.null(method)) {

    check.label(method, choices = available.deduplication.methods,
      labels = deduplication.labels, argname = "deduplication method",
      see = "dedup")

    if ((method == "cor") && (metadata$type != "continuous"))
      stop("method 'cor' may only be used with continuous data.")

    return(method)

  }#THEN
  else {

    return("cor")

  }#ELSE

}#CHECK.DEDUPLICATION.METHOD

check.deduplication.threshold = function(method, threshold) {

  if (method == "cor") {

    if (missing(threshold))
      threshold = 0.90
    else if (!is.probability(threshold))
      stop("the correlation threshold must be a number between 0 and 1.")

  }#THEN

  return(threshold)

}#CHECK.DEDUPLICATION.THRESHOLD
