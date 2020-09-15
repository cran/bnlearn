
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
    return(attr(strength, "threshold"))

  method = attr(strength, "method")

  if (method %in% c("test", "bootstrap", "bayes-factor")) {

    if (!is.probability(threshold))
      stop("the threshold must be a numeric value between 0 and 1.")

  }#THEN
  else if (method == "score"){

    if (!is.real.number(threshold) && !is.infinite(threshold))
      stop("the threshold must be numeric value.")

  }#THEN

  return(threshold)

}#CHECK.THRESHOLD

# check the cutpoints used to bin arc strengths.
check.cutpoints = function(cutpoints, strength, threshold, method) {

  if (!missing(cutpoints)) {

    if (method %in% c("test", "bootstrap", "bayes-factor")) {

      if (!is.probability(cutpoints))
        stop("the cutpoints must be numeric values between 0 and 1.")

    }#THEN
    else if (method == "score"){

      if (!is.real.number(cutpoints))
        stop("the cutpoints must be numeric values.")

    }#THEN

  }#THEN
  else {

    if (method == "test")
      cutpoints = unique(c(0, threshold/c(10, 5, 2, 1.5, 1), 1))
    else if (method %in% c("bootstrap", "bayes-factor"))
      cutpoints = unique(c(0, (1 - threshold)/c(10, 5, 2, 1.5, 1), 1))
    else if (method == "score") {

      # define a set of cut points using the quantiles from the empirical
      # distribution of the score deltas corresponding to significant arcs.
      significant = strength[strength < threshold]
      q = quantile(significant, 1 - c(0.50, 0.75, 0.90, 0.95, 1), names = FALSE)
      cutpoints = c(-Inf, threshold, unique(q), Inf)

    }#THEN

  }#ELSE

  return(sort(cutpoints))

}#CHECK.CUTPOINTS

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

