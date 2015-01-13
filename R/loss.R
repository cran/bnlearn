
loss.function = function(fitted, data, loss, extra.args, debug = FALSE) {

  if (loss %in% c("logl", "logl-g", "logl-cg")) {

    result = entropy.loss(fitted = fitted, data = data, debug = debug)

  }#THEN
  else if (loss == "pred") {

    result = classification.error(node = extra.args$target, fitted = fitted,
               prior = extra.args$prior, data = data, debug = debug)

  }#THEN
  else if (loss == "cor") {

    result = predictive.correlation(node = extra.args$target, fitted = fitted,
               data = data, debug = debug)

  }#THEN
  else if (loss == "mse") {

    result = mean.square.error(node = extra.args$target, fitted = fitted,
               data = data, debug = debug)

  }#THEN

  if (debug)
    cat("  @ total loss is", result$loss, ".\n")

  return(result)

}#LOSS.FUNCTION

# how to aggregate loss functions across the folds of cross-validation.
kfold.loss.postprocess = function(kcv, kcv.length, loss, extra.args, data) {

  if (loss == "cor") {

    if (all(is.na(unlist(sapply(kcv, "[", "loss"))))) {

      # this happens for root nodes, whose predictions are just the mean of
      # the response over the test set; propagate the NAs to the CV estimate.
      mean = NA

    }#THEN
    else {

      # match predicted and observed values.
      pred = unlist(lapply(kcv, "[[", "predicted"))
      obs = unlist(lapply(kcv, "[[", "observed"))
      # compute the predictive correlation.
      mean = cor(obs, pred)

    }#ELSE

  }#THEN
  else {

    # compute the mean of the observed values of the loss function, weighted
    # to account for unequal-length splits.
    mean = weighted.mean(unlist(sapply(kcv, '[', 'loss')), kcv.length)

  }#ELSE

  return(mean)

}#LOSS.POSTPROCESS

# predictive mean square error for gaussian networks.
mean.square.error = function(node, fitted, data, debug = FALSE) {

  pred = gaussian.prediction(node, fitted, data)

  return(list(loss = mean((data[, node] - pred)^2), predicted = pred,
    observed = data[, node]))

}#MEAN.SQUARE.ERROR

# predictive correlation for gaussian networks.
predictive.correlation = function(node, fitted, data, debug = FALSE) {

  pred = gaussian.prediction(node, fitted, data)

  if (length(fitted[[node]]$parents) == 0)
    return(list(loss = NA, predictions = pred, observed = data[, node]))
  else
    return(list(loss = cor(data[, node], pred), predicted = pred,
      observed = data[, node]))

}#PREDICTIVE.CORRELATION

# log-likelihood loss function for gaussian and discrete networks.
entropy.loss = function(fitted, data, keep = names(fitted), by.sample = FALSE,
    debug = FALSE) {

  list(loss = .Call("entropy_loss",
                    fitted = fitted,
                    data = data,
                    by.sample = by.sample,
                    keep = keep,
                    debug = debug))

}#ENTROPY.LOSS

# classification error as a loss function.
classification.error = function(node, fitted, prior = NULL, data, debug = FALSE) {

  if (is(fitted, c("bn.naive", "bn.tan")))
    pred = naive.classifier(node, fitted, prior, data)
  else
    pred = discrete.prediction(node, fitted, data)

  l = .Call("class_err",
            reference = minimal.data.frame.column(data, node),
            predicted = pred)

  if (debug)
    cat("  > classification error for node", node, "is", l, ".\n")

  return(list(loss = l, predicted = pred, observed = data[, node]))

}#CLASSIFICATION.ERROR

