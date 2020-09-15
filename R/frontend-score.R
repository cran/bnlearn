
# compute the score of a network.
network.score = function(x, data, type = NULL, ..., by.node = FALSE, debug = FALSE) {

  # check x's class.
  check.bn(x)
  # the original data set is needed.
  check.data(data)
  # check the network against the data.
  check.bn.vs.data(x, data)
  # check debug and by.node.
  check.logical(by.node)
  check.logical(debug)
  # no score if the graph is partially directed.
  if (is.pdag(x$arcs, names(x$nodes)))
    stop("the graph is only partially directed.")
  # check the score label.
  type = check.score(type, data)

  # expand and sanitize score-specific arguments.
  extra.args = check.score.args(score = type, network = x,
                 data = data, extra.args = list(...), learning = FALSE)
  # check that the score is decomposable when returning node contributions.
  if (by.node && !is.score.decomposable(type, extra.args))
    stop("the score is not decomposable, node terms are not defined.")

  # compute the node contributions to the network score.
  local = per.node.score(network = x, data = data, score = type,
            targets = names(x$nodes), extra.args = extra.args, debug = debug)

  if (by.node)
    return(local)
  else
    return(sum(local))

}#NETWORK.SCORE

# AIC method for class 'bn', an alias of score(..., type = "aic")
AIC.bn = function(object, data, ..., k = 1) {

  # check which type of data we are dealing with.
  type = data.type(data)

  # parameter sanitization done in the score() function.
  if (type %in% discrete.data.types)
    network.score(object, data = data, type = "aic", k = k, ...)
  else if (type == "continuous")
    network.score(object, data = data, type = "aic-g", k = k, ...)
  else if (type == "mixed-cg")
    network.score(object, data = data, type = "aic-cg", k = k, ...)

}#AIC.BN

# BIC method for class 'bn', an alias of score(..., type = "bic")
BIC.bn = function(object, data, ...) {

  # check which type of data we are dealing with.
  type = data.type(data)

  # parameter sanitization done in the score() function.
  if (type %in% discrete.data.types)
    network.score(object, data = data, type = "bic", ...)
  else if (type == "continuous")
    network.score(object, data = data, type = "bic-g", ...)
  else if (type == "mixed-cg")
    network.score(object, data = data, type = "bic-cg", ...)

}#BIC.BN

# logLik method for class 'bn', an alias of score(..., type = "loglik")
logLik.bn = function(object, data, ...) {

  # check which type of data we are dealing with.
  type = data.type(data)

  # parameter sanitization done in the score() function.
  if (type %in% discrete.data.types)
    network.score(x = object, data = data, type = "loglik", ...)
  else if (type == "continuous")
    network.score(x = object, data = data, type = "loglik-g", ...)
  else if (type == "mixed-cg")
    network.score(x = object, data = data, type = "loglik-cg", ...)

}#LOGLIK.BN

alpha.star = function(x, data, debug = FALSE) {

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

  alpha.star.backend(x = x, data = data, debug = debug)

}#ALPHA.STAR

# infer the direction of an ipothetic arc between two specified nodes.
choose.direction = function(x, arc, data, criterion = NULL, ..., debug = FALSE) {

  # check x's class.
  check.bn(x)
  # check the data are there.
  data.info = check.data(data)
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
  reset.test.counter()

  if (debug)
    cat("* testing", arc[1], "-", arc[2], "for direction.\n" )

  if (criterion %in% available.tests) {

    # sanitize the alpha threshold.
    alpha = check.alpha(list(...)$alpha, network = x)

    # sanitize B (the number of bootstrap/permutation samples).
    B = check.B(list(...)$B, criterion)

    # warn about unused arguments.
    check.unused.args(list(...), c("alpha", "B"))

    x = choose.direction.test(x, data = data, arc = arc,
          test = criterion, alpha = alpha, B = B, debug = debug,
          complete = data.info$complete.nodes)

  }#THEN
  else if (criterion %in% available.scores) {

    # expand and sanitize score-specific arguments.
    extra.args = check.score.args(score = criterion, network = x,
                   data = data, extra.args = list(...), learning = FALSE)

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

# compute the Bayes factor of two networks.
BF = function(num, den, data, score, ..., log = TRUE) {

  # check the two networks, individually and against each other.
  check.bn(num)
  check.bn(den)
  match.bn(num, den)
  nodes = names(num$nodes)
  # check the data.
  data.info = check.data(data)
  # check the networks against the data.
  check.bn.vs.data(num, data)
  check.bn.vs.data(den, data)
  # check the log argument.
  check.logical(log)
  # no score if at least one of the networks is partially directed.
  if (is.pdag(num$arcs, names(num$nodes)))
    stop("the graph in the numerator on the BF is only partially directed.")
  if (is.pdag(den$arcs, names(den$nodes)))
    stop("the graph in the denominator on the BF is only partially directed.")

  # make sure the score function is suitable for computing a Bayes factor.
  if (missing(score)) {

    if (data.info$type %in% discrete.data.types)
      score = "bde"
    else if (data.info$type %in% continuous.data.types)
      score = "bge"
    else if (data.info$type %in% mixed.data.types)
      score = "bic-cg"

  }#THEN
  else {

    score = check.score(score, data,
              allowed = c(available.discrete.bayesian.scores,
                          available.continuous.bayesian.scores,
                          grep("bic", available.scores, value = TRUE)))

  }#ELSE

  # expand and sanitize score-specific arguments.
  extra.args = check.score.args(score = score, network = num,
                 data = data, extra.args = list(...), learning = FALSE)

  # if a graph prior is used, this in not a Bayes factor any longer.
  if (!is.null(extra.args$prior) && extra.args$prior != "uniform")
    warning("using a non-uniform graph prior means this is not a Bayes factor.")

  # if the score is decomposable, compute the Bayes factor using only those
  # local distributions that differ between the two networks; otherwise
  # compute it on the whole network.
  if (is.score.decomposable(score, extra.args)) {

    different =
      sapply(nodes, function(n) {
        !setequal(num$nodes[[n]]$parents, den$nodes[[n]]$parents)
      })
    different = nodes[different]

  }#THEN
  else {

    different = nodes

  }#ELSE

  logBF.num = per.node.score(num, data = data, score = score,
                targets = different, extra.args = extra.args)
  logBF.den = per.node.score(den, data = data, score = score,
                targets = different, extra.args = extra.args)

  # compute the Bayes factor on the log-scale, and taking the difference between
  # local distributions before summing to minimise numeric problems.
  logBF = sum(logBF.num - logBF.den)

  return(ifelse(log, logBF, exp(logBF)))

}#BF
