
# predicted values for gaussian variables.
gaussian.prediction = function(node, fitted, data, debug = FALSE) {

  parents = fitted[[node]]$parents

  if (debug)
    cat("* predicting values for node ", node, ".\n", sep = "")

  if (length(parents) == 0) {

    .Call(call_gpred,
          fitted = fitted[[node]],
          ndata = nrow(data),
          debug = debug)

  }#THEN
  else {

    .Call(call_cgpred,
          fitted = fitted[[node]],
          parents = .data.frame.column(data, parents, drop = FALSE),
          debug = debug)

  }#ELSE

}#GAUSSIAN.PREDICTION

# predicted values for discrete networks.
discrete.prediction = function(node, fitted, data, prob = FALSE, debug = FALSE) {

  parents = fitted[[node]]$parents

  if (debug)
    cat("* predicting values for node ", node, ".\n", sep = "")

  if (length(parents) == 0) {

    .Call(call_dpred,
          fitted = fitted[[node]],
          ndata = nrow(data),
          prob = prob,
          debug = debug)

  }#THEN
  else {

    # if there is only one parent, get it easy.
    if (length(parents) == 1)
      config = .data.frame.column(data, parents)
    else
      config = configurations(.data.frame.column(data, parents), factor = FALSE)

    .Call(call_cdpred,
          fitted = fitted[[node]],
          parents = config,
          prob = prob,
          debug = debug)

  }#ELSE

}#DISCRETE.PREDICTION

# predicted values for conditional Gaussian networks.
mixedcg.prediction = function(node, fitted, data, debug = FALSE) {

  type = class(fitted[[node]])

  if (type == "bn.fit.dnode")
    discrete.prediction(node = node, fitted = fitted, data = data, debug = debug)
  else if (type == "bn.fit.gnode")
    gaussian.prediction(node = node, fitted = fitted, data = data, debug = debug)
  else {

    parents = fitted[[node]]$parents
    continuous.parents = names(which(sapply(parents, function(x) is.numeric(data[, x]))))
    discrete.parents = setdiff(parents, continuous.parents)

    # if there is only one parent, get it easy.
    if (length(discrete.parents) == 1)
      config = .data.frame.column(data, discrete.parents)
    else
      config = configurations(.data.frame.column(data, discrete.parents), factor = FALSE)

    .Call(call_ccgpred,
          fitted = fitted[[node]],
          configurations = config,
          parents = .data.frame.column(data, continuous.parents, drop = FALSE),
          debug = debug)

  }#ELSE

}#MIXEDCG.PREDICTION

# Naive Bayes and Tree-Augmented naive Bayes classifiers for discrete networks.
naive.classifier = function(training, fitted, prior, data, prob = FALSE,
    debug = FALSE) {

  # get the labels of the explanatory variables.
  nodes = names(fitted)
  # get the parents of each node, disregarding the training node.
  parents = sapply(fitted, function(x) {
    p = x$parents; return(p[p != training])
  })

  if (debug)
    cat("* predicting values for node ", training, ".\n", sep = "")

  .Call(call_naivepred,
        fitted = fitted,
        data = .data.frame.column(data, nodes, drop = FALSE),
        parents = match(parents, nodes),
        training = which(nodes == training),
        prior = prior,
        prob = prob,
        debug = debug)

}#NAIVE.CLASSIFIER

# maximum a posteriori predictions.
map.prediction = function(node, fitted, data, n, from, prob = FALSE,
    debug = FALSE) {

  .Call(call_mappred,
        node = node,
        fitted = fitted,
        data = data,
        n = as.integer(n),
        from = from,
        prob = prob,
        debug = debug)

}#MAP.PREDICTION

