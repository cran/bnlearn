
# check the list of networks passed to custom.strength().
check.customlist = function(custom, nodes) {

  objname = deparse(substitute(custom))

  # check that input is a list.
  if (!is(custom, "list"))
    stop(objname, " must be a list of objects of class 'bn', 'bn.fit' or of arc sets.")

  for (i in seq_along(custom)) {

    if (is(custom[[i]], c("bn", "bn.fit"))) {

      check.nodes(.nodes(custom[[i]]), graph = nodes,
        min.nodes = length(nodes), max.nodes = length(nodes))

    }#THEN
    else if (is(custom[[i]], "matrix")) {

      check.arcs(arcs = custom[[i]], nodes = nodes)

    }#THEN
    else {

      stop(objname, "[[", i, "]] is not an object of class 'bn', 'bn.fit' or an arc set.")

    }#ELSE

  }#FOR

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

# check a string that may be a score or a test label.
check.criterion = function(criterion, data) {

  if (!missing(criterion) && !is.null(criterion)) {

    # check and return errors from minimal.check.labels().
    check.label(criterion, choices = c(available.tests, available.scores),
      labels = c(test.labels, score.labels), argname = "criterion",
      see = "bnlearn-package")

  }#THEN
  else {

    # set the defaults using check.score() and check.test().
    if (criterion %in% available.tests)
      criterion = check.test(criterion, data)
    else if (criterion %in% available.scores)
      criterion = check.score(criterion, data)

  }#ELSE

  return(criterion)

}#CHECK.CRITERION

