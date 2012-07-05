
# compute the score of a network.
score = function(x, data, type = NULL, ..., debug = FALSE) {

  # check x's class.
  check.bn(x)
  # the original data set is needed.
  check.data(data)
  # check the network against the data.
  check.bn.vs.data(x, data)
  # check debug.
  check.logical(debug)
  # no score if the graph is partially directed.
  if (is.pdag(x$arcs, names(x$nodes)))
    stop("the graph is only partially directed.")
  # check the score label.
  type = check.score(type, data)

  # expand and sanitize score-specific arguments.
  extra.args = check.score.args(score = type, network = x,
                 data = data, extra.args = list(...))

  # compute the network score.
  network.score(network = x, data = data, score = type,
    extra.args = extra.args, debug = debug)

}#SCORE

# AIC method for class 'bn', an alias of score(..., type = "aic")
AIC.bn = function(object, data, ..., k = 1) {

  # parameter sanitization done in the score() function.
  if (is.data.discrete(data))
    score(object, data = data, type = "aic", k = k, ...)
  else
    score(object, data = data, type = "aic-g", k = k, ...)

}#AIC.BN

# logLik method for class 'bn', an alias of score(..., type = "loglik")
logLik.bn = function(object, data, ...) {

  # parameter sanitization done in the score() function.
  if (is.data.discrete(data))
    score(x = object, data = data, type = "loglik", ...)
  else
    score(x = object, data = data, type = "loglik-g", ...)

}#LOGLIK.BN

# infer the direction of an ipothetic arc between two specified nodes.
choose.direction = function(x, arc, data, criterion = NULL, ..., debug = FALSE) {

  # check x's class.
  check.bn(x)
  # check the data are there.
  check.data(data)
  # check the arc is there.
  check.nodes(nodes = arc, graph = x, min.nodes = 2, max.nodes = 2)
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
  else if (identical(criterion, "bootstrap")) {

    # nothing to do, move along.

  }#THEN
  else {

    criterion = check.criterion(criterion, data)

  }#ELSE

  # set the test/score counter.
  assign(".test.counter", 0, envir = .GlobalEnv)

  if (debug)
    cat("* testing", arc[1], "-", arc[2], "for direction.\n" )

  if (criterion %in% available.tests) {

    # sanitize the alpha threshold.
    alpha = check.alpha(list(...)$alpha, network = x)

    # sanitize B (the number of bootstrap/permutation samples).
    B = check.B(list(...)$B, criterion)

    # warn about unused arguments.
    check.unused.args(list(...), c("alpha", "B"))

    x = choose.direction.test(x, data = data, arc = arc, test = criterion,
          alpha = alpha, B = B, debug = debug)

  }#THEN
  else if (criterion %in% available.scores) {

    # expand and sanitize score-specific arguments.
    extra.args = check.score.args(score = criterion, network = x,
                   data = data, extra.args = list(...))

    x = choose.direction.score(x, data = data, arc = arc, score = criterion,
          extra.args = extra.args, debug = debug)

  }#ELSE
  else if (criterion == "bootstrap") {

    # expand and check bootstrap-specific arguments.
    extra.args = check.bootstrap.args(list(...), network = x, data = data)

    if (!is.null(extra.args$cpdag))
      check.logical(extra.args$cpdag)
    else
      extra.args$cpdag = TRUE

    x = choose.direction.boot(x, data = data, arc = arc,
          extra.args = extra.args, algorithm = extra.args[["algorithm"]],
          algorithm.args = extra.args[["algorithm.args"]],
          cpdag =  extra.args[["cpdag"]], debug = debug)

  }#THEN

  invisible(x)

}#CHOOSE.DIRECTION

