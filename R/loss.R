
loss.function = function(fitted, data, loss, extra.args, debug = FALSE) {

  if (loss == "logl") {

    result = entropy.loss(fitted = fitted, data = data)

  }#THEN
  else {

    result = predictive.loss(node = extra.args$target, fitted = fitted,
               prior = extra.args$prior, extra.args = extra.args,
               data = data, loss = loss)

  }#ELSE

  if (debug)
    cat("  @ effetive sample size is", result$effective.size, ".\n")

  return(result)

}#LOSS.FUNCTION

# compute the classification error.
clerr.loss = function(observed, predicted) {

  complete = !is.na(predicted) & !is.na(observed)
  how.many(predicted[complete] != observed[complete]) / how.many(complete)

}#CLERR.LOSS

# compute the (multiclass) F1 score.
f1.loss = function(observed, predicted) {

  # construct the confusion matrix...
  confusion.matrix = as.matrix(table(observed, predicted))
  # ... compute precision and recall...
  precision = diag(confusion.matrix) / colSums(confusion.matrix)
  recall = diag(confusion.matrix) / rowSums(confusion.matrix)
  # ... and combine them into the F1 or multiclass F1 score.
  f1 = 2 * precision * recall / (precision + recall)
  f1[is.na(f1)] = 0

  if (nlevels(observed) == 2)
    return(f1[1])
  else
    return(mean(f1))

}#F1.LOSS

# compute the area under the ROC curve.
auroc.loss = function(observed, probabilities) {

  auc = function(labels, probabilities) {

    # use the complete pairs like in the other losses.
    complete = !is.na(labels) & !is.na(probabilities)
    labels = labels[complete]
    probabilities = probabilities[complete]

    # identify the negative and positive classes...
    neg.class = levels(labels)[1]
    pos.class = levels(labels)[2]
    # ... the points at which FPR and TPR change...
    cutoffs = unique(c(sort(probabilities, decreasing = TRUE), 0))
    # ... and compute all the (FPR, TPR) points in the ROC curve.
    fast.pred = function(c) levels(labels)[(probabilities > c) + 1L]
    fp = sapply(cutoffs, function(c)
           sum((fast.pred(c) == pos.class) & (labels == neg.class)) )
    tp = sapply(cutoffs, function(c)
           sum((fast.pred(c) == pos.class) & (labels == pos.class)) )

    fpr = fp / sum(labels == neg.class)
    tpr = tp / sum(labels == pos.class)

    # keep those that are usable...
    complete = is.finite(fpr) & is.finite(tpr)
    fpr = fpr[complete]
    tpr = tpr[complete]
    # ... and give up if there are too few of them.
    if (sum(complete) < 2)
      return(NA)

    # compute the AUC.
    auc = 0
    for (i in 2:length(fpr))
      auc = auc + 0.5 * (fpr[i] - fpr[i - 1]) * (tpr[i] + tpr[i - 1])

    return(auc)

  }#AUC

  # if the prediction probabilties are not available, the loss is undefined.
  if (is.null(probabilities))
    return(NA)

  else if (nlevels(observed) == 2) {

    # return the AUROC of the positive class.
    auc(observed, probabilities[2, ])

  }#THEN
  else {

    # compute the one-vs-rest AUROCs...
    class.aucs = sapply(levels(observed), function(l) {

      adjusted.levels = levels(observed)
      adjusted.levels[adjusted.levels != l] = paste0("not", l)
      collapsed = factor(observed, labels = adjusted.levels)
      collapsed = relevel(collapsed, paste0("not", l))

      auc(collapsed, probabilities[l, ])

    })
    # ... and return the average.
    mean(class.aucs)

  }#ELSE

}#AUROC.LOSS

# compute the predictive correlation.
predcor.loss = function(observed, predicted) {

  complete = !is.na(predicted) & !is.na(observed)

  if ((cgsd(predicted[complete]) == 0) || (cgsd(observed[complete]) == 0))
    return(NA)
  else
    cor(observed, predicted, use = "complete.obs")

}#PREDCOR.LOSS

# compute the predictive mean square error.
mse.loss = function(observed, predicted) {

  mean((observed - predicted)^2, na.rm = TRUE)

}#MSE.LOSS

# how to aggregate loss functions across the folds of cross-validation.
kfold.loss.postprocess = function(kcv, kcv.length, loss, extra.args, data) {

  if (loss %in% c("pred", "pred-exact")) {

    clerr.loss(unlist(lapply(kcv, "[[", "observed")),
               unlist(lapply(kcv, "[[", "predicted")))

  }#THEN
  else if (loss == "f1") {

    f1.loss(unlist(lapply(kcv, "[[", "observed")),
            unlist(lapply(kcv, "[[", "predicted")))

  }#THEN
  else if (loss == "auroc") {

    auroc.loss(unlist(lapply(kcv, "[[", "observed")),
               do.call("cbind", lapply(kcv, "[[", "probabilities")))

  }#THEN
  else if (loss == "cor") {

    predcor.loss(unlist(lapply(kcv, "[[", "observed")),
                 unlist(lapply(kcv, "[[", "predicted")))

  }#THEN
  else if (loss == "mse") {

    mse.loss(unlist(lapply(kcv, "[[", "observed")),
             unlist(lapply(kcv, "[[", "predicted")))

  }#THEN
  else if (loss == "logl") {

    # compute the mean of the observed values of the loss function, weighted
    # to account for unequal-length splits and for missing values within splits.
    weighted.mean(x = unlist(sapply(kcv, '[', 'loss')),
                  w = unlist(sapply(kcv, '[', 'effective.size')))

  }#THEN

}#KFOLD.LOSS.POSTPROCESS

# how to aggregate the loss in hold-out cross-validation.
holdout.loss.postprocess = function(kcv, kcv.length, loss, extra.args, data) {

  if (loss %in% c("pred", "pred-exact")) {

    values = sapply(kcv, function(fold) {
      clerr.loss(fold$observed, fold$predicted)
    })

  }#THEN
  else if (loss == "f1") {

    values = sapply(kcv, function(fold) {
      f1.loss(fold$observed, fold$predicted)
    })

  }#THEN
  else if (loss == "auroc") {

    values = sapply(kcv, function(fold) {
      auroc.loss(fold$probabilities, fold$probabilities)
    })

  }#THEN
  else if (loss == "cor") {

    values = sapply(kcv, function(fold) {
      predcor.loss(fold$observed, fold$predicted)
    })

  }#THEN
  else if (loss == "mse") {

    values = sapply(kcv, function(fold) {
      mse.loss(fold$observed, fold$predicted)
    })

  }#THEN
  else if (loss == "logl") {

    values = sapply(kcv, '[[', 'loss')

  }#THEN

  if (loss == "logl") {

    # compute the mean of the observed values of the loss function, weighted
    # to account for unequal-length splits and for missing values within splits.
    mean = weighted.mean(x = values, w = sapply(kcv, '[[', 'effective.size'))
  }#THEN
  else {

    mean = mean(values)

  }#ELSE

  return(list(values = values, mean = mean))

}#HOLDOUT.LOSS.POSTPROCESS

# log-likelihood loss function.
entropy.loss = function(fitted, data) {

  attr(data, "metadata") = collect.metadata(data)
  l = loglikelihood(fitted, data = data)

  return(list(loss = - as.numeric(l) / attr(l, "nobs"),
              effective.size = attr(l, "nobs")))

}#ENTROPY.LOSS

# prepare predictions to compute the loss.
predictive.loss = function(node, fitted, prior = NULL, data, loss,
    extra.args) {

  # compute the predictions...
  if (loss == "pred-exact")
    pred = naive.classifier(node, fitted, prior, data, prob = TRUE)
  else
    pred = predict.backend(fitted = fitted, node = node, data = data,
             prob = TRUE, method = extra.args$predict,
             extra.args = extra.args$predict.args)
  # ... and extract the prediction probabilites, if available.
  prob = attr(pred, "prob")
  attr(pred, "prob") = NULL

  return(list(loss = NA, predicted = pred, observed = data[, node],
              probabilities = prob,
              effective.size = how.many(!is.na(pred) & !is.na(data[, node]))))

}#PREDICTIVE.LOSS

