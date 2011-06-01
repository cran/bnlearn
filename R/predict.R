
# predicted values for gaussian variables.
gaussian.prediction = function(node, fitted, data, debug = FALSE) {

  parents = fitted[[node]]$parents

  if (debug)
    cat("* predicting values for node ", node, ".\n", sep = "")

  if (length(parents) == 0) {

    .Call("gpred",
          fitted = fitted[[node]],
          data = nrow(data),
          debug = debug,
          PACKAGE = "bnlearn")

  }#THEN
  else {

    .Call("cgpred",
          fitted = fitted[[node]],
          data = minimal.data.frame.column(data, parents, drop = FALSE),
          debug = debug,
          PACKAGE = "bnlearn")

  }#ELSE

}#GAUSSIAN.PREDICTION

# predicted values for discrete networks.
discrete.prediction = function(node, fitted, data, debug = FALSE) {

  parents = fitted[[node]]$parents

  if (debug)
    cat("* predicting values for node ", node, ".\n", sep = "")

  if (length(parents) == 0) {

    .Call("dpred",
          fitted = fitted[[node]],
          data = minimal.data.frame.column(data, node),
          debug = debug,
          PACKAGE = "bnlearn")

  }#THEN
  else {

    # if there is only one parent, get it easy.
    if (length(parents) == 1)
      config = minimal.data.frame.column(data, parents)
    else
      config = configurations(minimal.data.frame.column(data, parents), factor = FALSE)

    .Call("cdpred",
          fitted = fitted[[node]],
          data = minimal.data.frame.column(data, node),
          parents = config,
          debug = debug,
          PACKAGE = "bnlearn")

  }#ELSE

}#DISCRETE.PREDICTION

# naive Bayes classifier for discrete networks.
naive.classifier = function(training, fitted, prior, data, debug = FALSE) {

  # get the labels of the explanatory variables.
  nodes = names(fitted)

  if (debug)
    cat("* predicting values for node ", training, ".\n", sep = "")

  .Call("naivepred",
        fitted = fitted,
        data = minimal.data.frame.column(data, nodes, drop = FALSE),
        training = which(nodes == training),
        prior = prior,
        debug = debug,
        PACKAGE = "bnlearn")

}#NAIVE.CLASSIFIER
