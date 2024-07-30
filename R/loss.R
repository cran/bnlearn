
loss.function = function(fitted, data, loss, extra.args, debug = FALSE) {

  if (loss %in% c("logl", "logl-g", "logl-cg")) {

    result = entropy.loss(fitted = fitted, data = data, debug = debug)

  }#THEN
  else if (loss %in% c("pred", "pred-exact", "pred-lw", "pred-lw-cg")) {

    result = classification.error(node = extra.args$target, fitted = fitted,
               prior = extra.args$prior, extra.args = extra.args,
               data = data, loss = loss, debug = debug)

  }#THEN
  else if (loss %in% c("cor", "cor-lw", "cor-lw-cg")) {

    result = predictive.correlation(node = extra.args$target, fitted = fitted,
               extra.args = extra.args, data = data, loss = loss, debug = debug)

  }#THEN
  else if (loss %in% c("mse", "mse-lw", "mse-lw-cg")) {

    result = mse.loss(node = extra.args$target, fitted = fitted,
               extra.args = extra.args, data = data, loss = loss, debug = debug)

  }#THEN

  if (debug)
    cat("  @ total loss is", result$loss, ".\n")

  return(result)

}#LOSS.FUNCTION

# how to aggregate loss functions across the folds of cross-validation.
kfold.loss.postprocess = function(kcv, kcv.length, loss, extra.args, data) {

  if (loss %in% c("cor", "cor-lw")) {

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
      mean = cor(obs, pred, use = "complete.obs")

    }#ELSE

  }#THEN
  else {

    # compute the mean of the observed values of the loss function, weighted
    # to account for unequal-length splits and for missing values within splits.
    mean = weighted.mean(x = unlist(sapply(kcv, '[', 'loss')),
                         w = unlist(sapply(kcv, '[', 'effective.size')))

  }#ELSE

  return(mean)

}#LOSS.POSTPROCESS

# predictive mean square error for gaussian networks.
mse.loss = function(node, fitted, data, loss, extra.args, debug = FALSE) {

  if (loss == "mse")
    pred = predict.backend(fitted = fitted, node = node, data = data,
             method = extra.args$predict, extra.args = extra.args$predict.args)
  else if (loss %in% c("mse-lw", "mse-lw-cg"))
    pred = predict.backend(fitted = fitted, node = node, data = data,
             method = "bayes-lw", extra.args = extra.args)

  effective.size = how.many(!is.na(pred) & !is.na(data[, node]))

  return(list(loss = mean((data[, node] - pred)^2, na.rm = TRUE),
              predicted = pred, observed = data[, node],
              effective.size = effective.size))

}#MSE.LOSS

# predictive correlation for gaussian networks.
predictive.correlation = function(node, fitted, extra.args, data, loss,
    debug = FALSE) {

  if (loss == "cor")
    pred = predict.backend(fitted = fitted, node = node, data = data,
             method = extra.args$predict, extra.args = extra.args$predict.args)
  else if (loss %in% c("cor-lw", "cor-lw-cg"))
    pred = predict.backend(fitted = fitted, node = node, data = data,
             method = "bayes-lw", extra.args = extra.args)

  effective.size = how.many(!is.na(pred) & !is.na(data[, node]))

  if (((loss == "cor") && length(fitted[[node]]$parents) == 0))
    return(list(loss = NA, predicted = pred, observed = data[, node]))
  else
    return(list(loss = cor(data[, node], pred, use = "complete.obs"),
                predicted = pred, observed = data[, node],
                effective.size = effective.size))

}#PREDICTIVE.CORRELATION

# log-likelihood loss function for gaussian and discrete networks.
entropy.loss = function(fitted, data, debug = FALSE) {

  attr(data, "metadata") = collect.metadata(data)
  l = loglikelihood(fitted, data = data, as.loss = TRUE, debug = debug)

  return(list(loss = - as.numeric(l) / attr(l, "nobs"),
              effective.size = attr(l, "nobs")))

}#ENTROPY.LOSS

# classification error as a loss function.
classification.error = function(node, fitted, prior = NULL, data, loss,
    extra.args, debug = FALSE) {

  if (loss == "pred-exact")
    pred = naive.classifier(node, fitted, prior, data)
  else if (loss == "pred")
    pred = predict.backend(fitted = fitted, node = node, data = data,
             method = extra.args$predict, extra.args = extra.args$predict.args)
  else if (loss %in% c("pred-lw", "pred-lw-cg"))
    pred = predict.backend(fitted = fitted, node = node, data = data,
             method = "bayes-lw", extra.args = extra.args)

  l = .Call(call_class_err,
            reference = .data.frame.column(data, node),
            predicted = pred)

  if (debug)
    cat("  > classification error for node", node, "is", l, ".\n")

  effective.size = how.many(!is.na(pred) & !is.na(data[, node]))

  return(list(loss = l, predicted = pred, observed = data[, node],
              effective.size = effective.size))

}#CLASSIFICATION.ERROR

