
# predicted values for gaussian variables.
gaussian.prediction = function(node, fitted, data) {

  parents = fitted[[node]]$parents

  if (length(parents) == 0) {

    .Call("gpred",
          fitted = fitted[[node]],
          data = nrow(data),
          PACKAGE = "bnlearn")

  }#THEN
  else {

    .Call("cgpred",
          fitted = fitted[[node]],
          data = minimal.data.frame.column(data, parents, drop = FALSE),
          PACKAGE = "bnlearn")

  }#ELSE

}#GAUSSIAN.PREDICTION

# predicted values for discrete networks.
discrete.prediction = function(node, fitted, data) {

  parents = fitted[[node]]$parents

  if (length(parents) == 0) {

    .Call("dpred",
          fitted = fitted[[node]],
          data = minimal.data.frame.column(data, node),
          PACKAGE = "bnlearn")

  }#THEN
  else {

    # if there is only one parent, get it easy.
    if (length(parents) == 1)
      config = minimal.data.frame.column(data, parents)
    else
      config = raw.configurations(minimal.data.frame.column(data, parents))

    .Call("cdpred",
          fitted = fitted[[node]],
          data = minimal.data.frame.column(data, node),
          parents = config,
          PACKAGE = "bnlearn")

  }#ELSE

}#DISCRETE.PREDICTION

# naive Bayes classifier for discrete networks.
naive.classifier = function(training, fitted, prior, data) {

  # get the prior distribution.
  levels.out = levels(data[, training])
  # get the labels of the explanatory variables.
  nodes = names(fitted)
  explanatory = nodes[nodes != training]

  pred = apply(data, 1, function(x) {

    prob = rep(1, length(prior))

    # multiply all the fitted probabilities.
    for (node in explanatory) 
      prob = prob * fitted[[node]]$prob[x[node], ]
    # add the prior distribution to the mix.
    prob = prob * prior
    # get which level is the most probable a posteriori.
    levels.out[which.max(prob)]

  })

  # convert the return value into a factor.
  pred = factor(pred, levels.out)

  return(pred)

}#NAIVE.CLASSIFIER
