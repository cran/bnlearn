
# measure the strength of the arcs in a directed graph.
arc.strength = function(x, data, criterion = NULL, ..., debug = FALSE) {

  # check x's class.
  check.bn(x)
  # arc strength is undefined in partially directed graphs.
  if (is.pdag(x$arcs, names(x$nodes)) && !identical(criterion, "bootstrap"))
    stop("the graph is only partially directed.")
  # check the data are there.
  check.data(data)
  # check the network against the data
  check.bn.vs.data(x, data)
  # check debug.
  check.logical(debug)
  # check criterion.
  if (is.null(criterion)) {

    # if no criterion is specified use either the default one or the
    # one used by the learning algorithm.
    if (x$learning$test == "none")
      criterion = check.test(criterion, data)
    else
      criterion = x$learning$test

  }#THEN
  else if (!identical(criterion, "bootstrap")) {

    criterion = check.criterion(criterion, data)

  }#THEN

  # expand and sanitize score-specific arguments and the alpha threshold.
  if (criterion %in% available.tests) {

    # sanitize the alpha threshold.
    alpha = check.alpha(list(...)$alpha, network = x)

    # sanitize B (the number of bootstrap/permutation samples).
    B = check.B(list(...)$B, criterion)

    # warn about unused arguments.
    check.unused.args(list(...), c("alpha", "B"))

    res = arc.strength.test(network = x, data = data, alpha = alpha,
            test = criterion, B = B, debug = debug)

    # add extra information for strength.plot().
    res = structure(res, mode = "test", threshold = alpha)

  }#THEN
  else if (criterion %in% available.scores) {

    # expand and sanitize score-specific arguments.
    extra.args = check.score.args(score = criterion, network = x,
                   data = data, extra.args = list(...))

    res = arc.strength.score(network = x, data = data, score = criterion,
            extra = extra.args, debug = debug)

    # add extra information for strength.plot().
    res = structure(res, mode = "score", threshold = 0)

  }#THEN
  else if (criterion == "bootstrap") {

    # expand and check bootstrap-specific arguments.
    extra.args = check.bootstrap.args(list(...), network = x, data = data)

    res = arc.strength.boot(data = data, R = extra.args$R, 
            m = extra.args$m, algorithm = extra.args[["algorithm"]], 
            algorithm.args = extra.args[["algorithm.args"]], arcs = x$arcs, 
            debug = debug)

    # add extra information for strength.plot(), and drop the column
    # with the direction confidence.
    res = structure(res[, 1:3], mode = "bootstrap", threshold = 0.5)

  }#THEN

  return(structure(res, class = c("bn.strength", class(res))))

}#ARC.STRENGTH

