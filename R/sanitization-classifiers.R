
# sanitize the extra arguments passed to Bayesian classifiers.
check.classifier.args = function(method, data, training, explanatory,
    extra.args) {

  # check the label of the mutual information estimator.
  if (has.argument(method, "estimator", learning.extra.args))
    extra.args[["estimator"]] =
      check.mi.estimator(extra.args[["estimator"]], data = data)

  # check the node to use the root of the tree (if not specified pick the first
  # explanatory variable assuming natural ordering).
  if (has.argument(method, "root", learning.extra.args)) {

    if (!is.null(extra.args[["root"]]))
      check.nodes(extra.args[["root"]], graph = explanatory, max.nodes = 1)
    else
      extra.args[["root"]] = explanatory[1]

  }#THEN

  # warn about and remove unused arguments.
  extra.args = check.unused.args(extra.args, learning.extra.args[[method]])

  return(extra.args)

}#CHECK.CLASSIFIER.ARGS

# check a prior distribution against the observed variable.
check.classifier.prior = function(prior, training) {

  if (missing(prior) || is.null(prior)) {

    # use the empirical probabilities in the fitted network, or a flat prior
    # as a last resort.
    if (is(training, c("bn.fit.dnode", "bn.fit.onode")))
      prior = training$prob
    else
      prior = rep(1, nlevels(training))

  }#THEN
  else {

    if (is(training, c("bn.fit.dnode", "bn.fit.onode")))
      nlvls = dim(training$prob)[1]
    else
      nlvls = nlevels(training)

    if (length(prior) != nlvls)
      stop("the prior distribution and the training variable have a different number of levels.")
    if (!is.nonnegative.vector(prior))
      stop("the prior distribution must be expressed as a probability vector.")

    # make sure the prior probabilities sum to one.
    prior = prior / sum(prior)

  }#ELSE

  return(prior)

}#CHECK.CLASSIFIER.PRIOR

