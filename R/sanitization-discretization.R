
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

  if (method == "hartemink") {

    if (!is.null(extra.args$idisc)) {

      other.methods = available.discretization.methods[available.discretization.methods != "hartemink"]

      check.label(extra.args$idisc, choices = other.methods,
        labels = discretization.labels, argname = "initial discretization method",
        see = "discretize")

    }#THEN
    else {

      # default to quantile discretization as per Hartemink's recommendation.
      extra.args$idisc = "quantile"

    }#ELSE

    if (!is.null(extra.args$ibreaks)) {

      if (length(extra.args$ibreaks) == 1)
        extra.args$ibreaks = rep(extra.args$ibreaks, ncol(data))
      else if (length(extra.args$ibreaks) != ncol(data))
        stop("the 'ibreaks' vector must have an element for each variable in the data.")
      if (!is.positive.vector(extra.args$ibreaks))
        stop("the number of initial breaks must be a (vector of) positive integer number(s).")
      if (any(extra.args$ibreaks < breaks))
        stop("insufficient number of initial breaks, need at least ", breaks + 1, ".")
      else if (any(extra.args$ibreaks == breaks))
        warning("the initial number of breaks is identical to the final number of breaks.")

    }#THEN
    else {

      ndata = nrow(data)

      if (ndata > 500)
        extra.args$ibreaks = 50
      else if (ndata > 100)
        extra.args$ibreaks = 20
      else if (ndata > 50)
        extra.args$ibreaks = 10
      else if (ndata > 10)
        extra.args$ibreaks = 5
      else
        extra.args$ibreaks = ndata

      extra.args$ibreaks = rep(extra.args$ibreaks, ncol(data))

    }#ELSE

  }#THEN

  # warn about unused arguments.
  check.unused.args(extra.args, discretization.extra.args[[method]])

  return(extra.args)

}#CHECK.DISCRETIZATION.ARGS

