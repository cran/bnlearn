
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
