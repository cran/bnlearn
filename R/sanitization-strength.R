
# check the list of networks passed to custom.strength().
check.customlist = function(custom, nodes) {

  # check
  if (!is(custom, "list"))
    stop("networks must be a list of objects of class 'bn' or of arc sets.")
  if (!all(sapply(custom, function(x) { is(x, "bn") || is(x, "matrix") })))
    stop("x must be a list of objects of class 'bn' or of arc sets.")

  validate = function(custom, nodes) {

    if (is(custom, "bn")) {

      check.nodes(names(custom$nodes), graph = nodes, min.nodes = length(nodes),
        max.nodes = length(nodes))

    }
    else if (is(custom, "matrix")) {

      check.arcs(arcs = custom, nodes = nodes)

    }#THEN
    else {

      stop("x must be a list of objects of class 'bn' or of arc sets.")

    }#ELSE

  return(TRUE)

  }#VALIDATE

  if (!all(sapply(custom, validate, nodes = nodes)))
    stop("x must be a list of objects of class 'bn' or of arc sets.")

}#CHECK.CUSTOMLIST

# check an object of class bn.strength.
check.bn.strength = function(strength, nodes, valid = available.strength.methods) {

  # check that the object is there and with the right class.
  if (missing(strength))
    stop("an object of class 'bn.strength' is required.")
  if (!is(strength, "bn.strength"))
    stop(sprintf("%s must be an object of class 'bn.strength'.",
           deparse(substitute(strength))))
  # check the object structure.
  if (ncol(strength) %!in% 3:4)
    stop("objects of class 'bn.strength' must have 3 or 4 columns.")
  if (!identical(names(strength), c("from", "to", "strength")) &&
      !identical(names(strength), c("from", "to", "strength", "direction")))
    stop("objects of class 'bn.strength' must be data frames with column names ",
         "'from', 'to', 'strength' and (optionally) 'direction'.")
  if (any(c("method", "threshold") %!in% names(attributes(strength))))
    stop("objects of class 'bn.strength' must have a 'method' and a 'strength' attributes.")
  # check the estimation method.
  if (attr(strength, "method") %!in% valid)
    check.label(attr(strength, "method"), choices = valid,
      argname = "strength estimation methods")
  # check the consistency with a network's node set.
  if (!missing(nodes))
    check.arcs(strength[, c("from", "to"), drop = FALSE], nodes)

}#CHECK.BN.STRENGTH

# sanitize the threshold value.
check.threshold = function(threshold, strength) {

  if (missing(threshold))
    threshold = attr(strength, "threshold")
  else {

    s = strength[, "strength"]

    if (!is.numeric(threshold) || (length(threshold) != 1) || is.nan(threshold))
      stop("the threshold must be a numeric value.")
    if ((threshold < min(s)) || (threshold > max(s)))
      warning("the threshold is outside the range of the strength values.")

  }#ELSE

  return(threshold)

}#CHECK.THRESHOLD

