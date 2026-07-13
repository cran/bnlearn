
# prediction imputation backend.
predict.backend = function(fitted, node, data, cluster = NULL, method,
    extra.args, prob = FALSE, debug = FALSE) {

  # individual incomplete observations are imputed independently of each other,
  # so they can be processed in parallel.
  if (!is.null(cluster))
    splits = parallel::clusterSplit(cluster, seq(nrow(data)))
  else
    splits = list(seq(nrow(data)))

  # choose the right function for the prediction method.
  if (method == "parents") {

    fun = parents.prediction

  }#THEN
  else if (method == "bayes-lw") {

    fun = likelihood.weighting.prediction

  }#THEN
  else if (method == "exact") {

    if (is(fitted, c("bn.fit.dnet", "bn.fit.onet", "bn.fit.donet"))) {

      if (is(fitted, c("bn.fit.onet", "bn.fit.donet")))
        warning("the gRain package does not support ordinal networks, disregarding the ordering of the levels.")

      fun = exact.discrete.prediction

    }#THEN
    else if (is(fitted, c("bn.fit.gnet"))) {

      fun = exact.gaussian.prediction

    }#THEN
    else
      stop("'bn.fit.cgnet' networks are not supported.")

  }#THEN

  # compute the predicted values.
  values = smartSapply(cluster, splits, function(ids) {
    fun(node = node, fitted = fitted, data = data[ids, , drop = FALSE],
        extra.args = extra.args, prob = prob, debug = debug)
  })

  if (!is.null(cluster)) {

    predicted = do.call("c", values)

    # reassemble and attach the classification probabilities.
    if (prob) {

      clprobs = lapply(values, attr, "prob")
      attr(predicted, "prob") = do.call(cbind, clprobs)

    }#THEN

  }#THEN
  else {

    predicted = values[[1]]

  }#THEN

  return(predicted)

}#PREDICT.BACKEND

# predicted values from the parents of a node; the underlying C code handles
# discrete, Gaussian, conditional Gaussian and zero-inflated nodes.
parents.prediction = function(node, fitted, data, extra.args, prob = FALSE,
    debug = FALSE) {

  if (debug)
    cat("* predicting values for node ", node, ".\n", sep = "")

  .Call(call_predict_from_parents,
        fitted = fitted,
        node = node,
        data = data,
        prob = prob,
        debug = debug)

}#PARENTS.PREDICTION

# Naive Bayes and Tree-Augmented naive Bayes classifiers for discrete networks.
naive.classifier = function(training, fitted, prior, data, prob = FALSE,
    debug = FALSE) {

  if (debug)
    cat("* predicting values for node ", training, ".\n", sep = "")

  .Call(call_naivepred,
        fitted = fitted,
        data = data,
        training = which(names(fitted) == training),
        prior = prior,
        prob = prob,
        debug = debug)

}#NAIVE.CLASSIFIER

# predict a node from an incomplete set of predictors by imputing it while
# averaging over the unobserved predictors, using the given imputation backend;
# the backend also returns the prediction probabilities, if requested.
predict.by.imputing = function(node, fitted, data, extra.args, impute.fun,
    prob = FALSE, debug = FALSE) {

  # if the predictors include the target variable, erase its values; otherwise
  # introduce it as a column of NAs (a factor for discrete targets) to be imputed.
  if (node %in% names(data))
    data[, node][] = NA
  else {

    data[, node] = rep(NA_real_, nrow(data))
    if (is(fitted[[node]], c("bn.fit.dnode", "bn.fit.onode")))
      data[, node] = factor(data[, node],
                            levels = dimnames(fitted[[node]]$prob)[[1]])

  }#ELSE

  imputed = impute.fun(fitted = fitted, data = data, extra.args = extra.args,
              restrict.from = extra.args$from, restrict.to = node, prob = prob,
              debug = debug)

  # subsetting drops attributes: carry over the posterior probabilities.
  predicted = imputed[, node]
  if (prob)
    attr(predicted, "prob") = attr(imputed, "prob")

  return(predicted)

}#PREDICT.BY.IMPUTING

# maximum a posteriori predictions.
likelihood.weighting.prediction = function(node, fitted, data, extra.args,
    prob = FALSE, debug = FALSE) {

  complete.nodes = attr(data, "metadata")$complete.nodes
  complete.predictors = complete.nodes[extra.args$from]

  # if the data are incomplete, the set of predictors may differ between
  # observations, so predict by imputing the target instead.
  if (!all(complete.predictors))
    return(predict.by.imputing(node = node, fitted = fitted, data = data,
             extra.args = extra.args,
             impute.fun = impute.backend.likelihood.weighting, prob = prob,
             debug = debug))

  .Call(call_predict_map_lw,
        node = node,
        fitted = fitted,
        data = data,
        n = as.integer(extra.args$n),
        from = extra.args$from,
        prob = prob,
        debug = debug)

}#LIKELIHOOD.WEIGHTING.PREDICTION

# prediction using exact inference in discrete networks.
exact.discrete.prediction = function(node, fitted, data, extra.args,
    restrict.from = NULL, restrict.to = NULL, prob = FALSE, debug = FALSE) {

  complete.nodes = attr(data, "metadata")$complete.nodes
  complete.predictors = complete.nodes[extra.args$from]

  # if the data are incomplete, the set of predictors may differ between
  # observations, so predict by imputing the target instead.
  if (!all(complete.predictors))
    return(predict.by.imputing(node = node, fitted = fitted, data = data,
             extra.args = extra.args, impute.fun = impute.backend.exact,
             prob = prob, debug = debug))

  # check whether gRain is loaded.
  check.and.load.package("gRain")

  jtree = from.bn.fit.to.grain(fitted, compile = FALSE)
  compact.jtree = gRbase::compile(jtree, propagate = FALSE)

  # check the size of the conditional probability table implied by the query.
  node.nlevels = sapply(fitted, function(x) dim(x$prob)[1])
  from = extra.args$from
  query.size = prod(node.nlevels[c(node, from)])

  # if the conditional probability table is small (20MB), check whether the
  # junction tree with all the query nodes in the root clique (possibly more
  # memory, faster inference) has similar size to the most compact one (less
  # memory, slower inference).
  if (query.size <= 20 * 1024^2 / 8) {

    targeted.jtree =
      gRbase::compile(jtree, root = c(node, from), propagate = FALSE)

    # if both the conditional probability table and the targeted junction tree
    # are small enough, materialize the conditional probability table and reuse
    # the code that predicts from the parents; otherwise, use the compact tree
    # and iterate over the observations.
    if ((tree.size(targeted.jtree, node.nlevels) <=
          2 * tree.size(compact.jtree, node.nlevels))) {

      predval = targeted.exact.discrete.prediction(jtree = targeted.jtree,
                  node = node, data = data, from = from, prob = prob,
                  debug = debug)
    }#ELSE
    else {

      predval = compact.exact.discrete.prediction(jtree = compact.jtree,
                  node = node, data = data, from = from, prob = prob,
                  debug = debug)

    }#ELSE

  }#THEN
  else {

    # if the conditional probability table is too large, use the compact
    # tree iterate over the observations.
    predval = compact.exact.discrete.prediction(jtree = compact.jtree,
                node = node, data = data, from = from, prob = prob,
                debug = debug)

  }#ELSE

  return(predval)

}#EXACT.DISCRETE.PREDICTION

# prediction using exact inference to create a conditional probability table
# to predict from.
targeted.exact.discrete.prediction = function(jtree, node, data, from,
    prob = FALSE, debug = FALSE) {

  # generate the conditional probability table of the target variable given
  # the predictors...
  cpt = grain.query(jtree, nodes = c(node, from), type = "conditional")
  # ... embed it in a minimal bn.fit object...
  fitted.node = structure(list(node = node, parents = from, prob = cpt),
                          class = "bn.fit.dnode")
  fitted.network = structure(list(fitted.node), names = node)
  # ... and perform prediction treating predictors as parents.
  predval = parents.prediction(node = node, fitted = fitted.network,
              data = data, prob = prob, debug = debug)

  return(predval)

}#TARGETED.EXACT.DISCRETE.PREDICTION

# prediction using exact inference in discrete networks, iterating over
# observations and propagating the predictors as evidence.
compact.exact.discrete.prediction = function(jtree, node, data, from,
    prob = FALSE, debug = FALSE) {

  # the predicted values are a factor with the target node's levels, taken
  # from the junction tree since the target node is not necessarily available.
  predval = factor(rep(NA, nrow(data)), levels = jtree$universe$levels[[node]])
  # the prediction probabilities, if we are to return them.
  if (prob)
    classprob = matrix(NA, nrow = nlevels(predval), ncol = nrow(data),
                  dimnames = list(levels(predval), NULL))

  if (debug)
    cat("* predicting values for node ", node, ".\n", sep = "")

  # for each observation...
  for (i in seq(nrow(data))) {

    # ... set the predictors as evidence...
    jpred = gRain::setEvidence(jtree, nodes = from,
              states = sapply(data[i, from, drop = FALSE], as.character))
    # ... check that it is possible to observe the evidence...
    if (gRain::pEvidence(jpred) <= sqrt(.Machine$double.eps)) {

      predval[i] = NA
      if (prob)
        classprob[, i] = NA
      next

    }#THEN
    # ... get the probability distribution of the node being predicted...
    predprob = grain.query(jpred, nodes = node, type = "marginal")[[node]]
    # ... and choose the level with the highest probability, breaking ties
    # randomly.
    all.maxima = names(predprob[predprob == max(predprob)])

    if (length(all.maxima) == 1)
      predval[i] = all.maxima
    else
      predval[i] = sample(all.maxima, size = 1)

    # also save the prediction probabilities.
    if (prob)
      classprob[, i] = predprob

    if (debug) {

      cat("  > prediction for observation ", i, " is '",
          as.character(predval[i]), "' with probabilities:\n", sep = "")
      cat("  ", predprob, "\n", sep = "  ")

    }#THEN

  }#FOR

  # attach the prediction probabilities to the return value if requested to.
  if (prob)
    attr(predval, "prob") = classprob

  return(predval)

}#COMPACT.EXACT.DISCRETE.PREDICTION

# prediction using exact inference in gaussian networks.
exact.gaussian.prediction = function(node, fitted, data, extra.args,
    prob = FALSE, debug = FALSE) {

  complete.nodes = attr(data, "metadata")$complete.nodes
  complete.predictors = complete.nodes[extra.args$from]

  # if the data are incomplete, the set of predictors may differ between
  # observations, so predict by imputing the target instead.
  if (!all(complete.predictors))
    return(predict.by.imputing(node = node, fitted = fitted, data = data,
             extra.args = extra.args, impute.fun = impute.backend.exact,
             prob = prob, debug = debug))

  # get the global distribution...
  mvn = gbn2mvnorm.backend(fitted)
  # ... reduce it to the target variable and the predictors...
  from = extra.args$from
  mu = mvn$mu[c(node, from), drop = FALSE]
  sigma = mvn$sigma[c(node, from), c(node, from), drop = FALSE]
  # ... check that the model is not singular...
  if (anyNA(mu) || anyNA(sigma))
    return(rep(NA_real_, nrow(data)))
  # ... create the reduced network that will produce the predictions...
  if (length(from) == 0) {

    prediction.network = model2network(paste0("[", node, "]"))

  }#THEN
  else {

    prediction.network = model2network(paste0(
      paste0("[", node, "|", paste0(from, collapse = ":"), "]"),
      paste0("[", from, "]", collapse = "")
    ))

  }#ELSE
  distribution = mvnorm2gbn.backend(prediction.network, mu = mu, sigma = sigma)
  # ... and predict.
  parents.prediction(node = node, fitted = distribution, data = data,
    debug = debug)

}#EXACT.GAUSSIAN.PREDICTION
