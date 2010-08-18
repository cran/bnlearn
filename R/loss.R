
loss.function = function(fitted, data, loss, extra.args, debug = FALSE) {

  result = 0
  nodes = names(fitted)

  if (loss == "logl") {

    result = discrete.loss(nodes = nodes, fitted = fitted, data = data,
               debug = debug)

  }#THEN
  else if (loss == "logl-g") {

    result = gaussian.loss(nodes = nodes, fitted = fitted, data = data,
               debug = debug)

  }#THEN
  else if (loss == "pred") {

    result = classification.error(node = extra.args$target, fitted = fitted,
      data = data, debug = debug)

  }#THEN

  if (debug)
    cat("  @ total loss is", result, ".\n")

  return(result)

}#LOSS.FUNCTION

# log-likelihood loss function for gaussian networks.
gaussian.loss = function(nodes, fitted, data, debug = FALSE) {

  sum(sapply(nodes, gaussian.loss.node, fitted = fitted, data = data,
    debug = debug))

}#GAUSSIAN.LOSS

# log-likelihood loss function for nodes in a gaussian network.
gaussian.loss.node = function(node, fitted, data, debug = FALSE) {

  parents = fitted[[node]]$parents

  if (length(parents) == 0) {

    l = .Call("gloss",
              fitted = fitted[[node]],
              data = minimal.data.frame.column(data, node),
              PACKAGE = "bnlearn")

  }#THEN
  else {

    l = .Call("cgloss",
              fitted = fitted[[node]],
              data = minimal.data.frame.column(data, c(node, parents)),
              PACKAGE = "bnlearn")

  }#ELSE

  if (debug)
    cat("  > log-likelihood loss for node", node, "is", l, ".\n")

  return(l)

}#GAUSSIAN.LOSS

# log-likelihood loss function for discrete networks.
discrete.loss = function(nodes, fitted, data, debug = FALSE) {

  sum(sapply(nodes, discrete.loss.node, fitted = fitted, data = data,
    debug = debug))

}#DISCRETE.LOSS

# log-likelihood loss function for nodes in a discrete network.
discrete.loss.node = function(node, fitted, data, debug = FALSE) {

  parents = fitted[[node]]$parents

  if (length(parents) == 0) {

    l = .Call("dloss",
              fitted = fitted[[node]],
              data = minimal.data.frame.column(data, node),
              node = node,
              PACKAGE = "bnlearn")

  }#THEN
  else {

    # if there is only one parent, get it easy.
    if (length(parents) == 1)
      config = minimal.data.frame.column(data, parents)
    else
      config = raw.configurations(minimal.data.frame.column(data, parents))

    l = .Call("cdloss",
              fitted = fitted[[node]],
              data = minimal.data.frame.column(data, node),
              parents = config,
              node = node,
              PACKAGE = "bnlearn")

  }#ELSE

  if (debug)
    cat("  > log-likelihood loss for node", node, "is", l, ".\n")

  return(l)

}#DISCRETE.LOSS

# classification error as a loss function.
classification.error = function(node, fitted, data, debug = FALSE) {

  l = .Call("class_err",
            reference = minimal.data.frame.column(data, node),
            predicted = discrete.prediction(node, fitted, data),
            PACKAGE = "bnlearn")

  if (debug)
    cat("  > classification error for node", node, "is", l, ".\n")

  return(l)

}#CLASSIFICATION.ERROR