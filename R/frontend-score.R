
# compute the score of a network.
network.score = function(x, data, type = NULL, ..., by.node = FALSE, debug = FALSE) {

  # check x's class.
  check.bn(x)
  # the original data set is needed.
  data = check.data(data, allow.missing = TRUE)
  # check the network against the data.
  check.bn.vs.data(x, data)
  # check debug and by.node.
  check.logical(by.node)
  check.logical(debug)
  # no score if the graph is partially directed.
  if (is.pdag(x$arcs, names(x$nodes)))
    stop("the graph is only partially directed.")
  # check the score label.
  type = check.score(type, data = data)
  # check whether the network is valid for the method.
  check.arcs.against.assumptions(x$arcs, data, type)

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

  # argument sanitization is performed in the score() function.
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

  # argument sanitization is performed in the score() function.
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

  # argument sanitization is performed in the score() function.
  if (type %in% discrete.data.types)
    network.score(x = object, data = data, type = "loglik", ...)
  else if (type == "continuous")
    network.score(x = object, data = data, type = "loglik-g", ...)
  else if (type == "mixed-cg")
    network.score(x = object, data = data, type = "loglik-cg", ...)

}#LOGLIK.BN

# estimate a reasonable guess at the best imaginary sample size.
alpha.star = function(x, data, debug = FALSE) {

  # check x's class.
  check.bn(x)
  # the original data set is needed.
  check.data(data, allowed.types = discrete.data.types)
  # check the network against the data.
  check.bn.vs.data(x, data)
  # check debug.
  check.logical(debug)
  # no score if the graph is partially directed.
  if (is.pdag(x$arcs, names(x$nodes)))
    stop("the graph is only partially directed.")

  alpha.star.backend(x = x, data = data, debug = debug)

}#ALPHA.STAR

# compute the Bayes factor of two networks.
BF = function(num, den, data, score, ..., log = TRUE) {

  # check the two networks, individually and against each other.
  check.bn(num)
  check.bn(den)
  match.bn(num, den)
  nodes = names(num$nodes)
  # check the data.
  data = check.data(data, allow.missing = TRUE)
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
  data.type = attr(data, "metadata")$type

  if (missing(score)) {

    if (data.type %in% discrete.data.types)
      score = "bde"
    else if (data.type %in% continuous.data.types)
      score = "bge"
    else if (data.type %in% mixed.data.types)
      score = "bic-cg"

  }#THEN
  else {

    score = check.score(score, data = data,
              allowed = c(available.discrete.bayesian.scores,
                          available.continuous.bayesian.scores,
                          available.omnibus.scores,
                          grep("bic", available.scores, value = TRUE)))

  }#ELSE

  # check whether the networks are valid for the score.
  check.arcs.against.assumptions(num$arcs, data, score)
  check.arcs.against.assumptions(den$arcs, data, score)

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
