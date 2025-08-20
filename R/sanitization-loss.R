
# check loss functions' labels.
check.loss = function(loss, data, bn) {

  if (!is.null(loss)) {

    # check the loss function.
    check.label(loss, choices = loss.functions, labels = loss.labels,
      argname = "loss function", see = "bn.cv")

    if (loss %in% classifiers.loss.functions) {

      if ((is.character(bn) && (bn %!in% classification.algorithms)) ||
          (!is.character(bn) && !is(bn, available.classifiers)))
        stop("loss function '", loss, "' may only be used with classifiers.")

    }#THEN

    return(loss)

  }#THEN
  else {

    if ((is.character(bn) && (bn %in% classification.algorithms)) ||
         is(bn, available.classifiers))
      return("pred-exact")

    return("logl")

  }#ELSE

}#CHECK.LOSS

# sanitize the extra arguments passed to loss functions.
check.loss.args = function(loss, bn, nodes, data, extra.args) {

  # target node that that the loss is computed for.
  if (has.argument(loss, "target", loss.extra.args))
    extra.args[["target"]] =
      check.loss.target(target = extra.args[["target"]], loss = loss,
        data = data, bn = bn, nodes = nodes)

  # method used for prediction and its optional arguments.
  if (has.argument(loss, "predict", loss.extra.args))
    extra.args[["predict"]] =
      check.prediction.method(method = extra.args[["predict"]])
  if (has.argument(loss, "predict.args", loss.extra.args)) {

    dummy.bn = structure(vector(ncol(data), mode = "list"),
                 names = colnames(data), class = "bn.fit")

    extra.args[["predict.args"]] =
      check.prediction.extra.args(method = extra.args[["predict"]],
        extra.args = extra.args[["predict.args"]],
        node = extra.args[["target"]], fitted = dummy.bn, data = data)

  }#THEN

  # prior distribution for the target variable of a classifier.
  if (has.argument(loss, "prior", loss.extra.args))
      extra.args[["prior"]] =
        check.classifier.prior(prior = extra.args[["prior"]],
          training = data[, extra.args[["target"]]])

  # number of particles for Monte Carlo-based losses.
  if (has.argument(loss, "n", loss.extra.args)) {

    if (!is.null(extra.args[["n"]])) {

      if (!is.positive.integer(extra.args[["n"]]))
        stop("the number of observations to be sampled must be a positive integer number.")

    }#THEN
    else {

      extra.args[["n"]] = 500

    }#ELSE

  }#THEN

  # which nodes to predict from (all of them unless told otherwise).
  if (has.argument(loss, "from", loss.extra.args)) {

    if (!is.null(extra.args[["from"]]))
      check.nodes(extra.args[["from"]], graph = names(data), min.nodes = 1)
    else
      extra.args[["from"]] = setdiff(names(data), extra.args[["target"]])

  }#THEN

  # warn about and remove unused arguments.
  extra.args = check.unused.args(extra.args, loss.extra.args[[loss]])

  return(extra.args)

}#CHECK.LOSS.ARGS

# check the target node the loss is computed for.
check.loss.target = function(target, loss, data, bn, nodes) {

  if (!is.null(target)) {

    if (!is.string(target) || (target %!in% nodes))
      stop("target node must be a single, valid node label for the network.")

    # check the target has the right data type.
    if (loss %in% c("cor", "mse"))
      if (!is(data[, target], "numeric"))
        stop("the target node must be a continuous variable.")
    if (loss %in% c("pred", "pred-exact", "f1"))
      if (!is(data[, target], "factor"))
        stop("the target node must be a factor.")

  }#THEN
  else {

    # the target node is obvious for classifiers.
    if (is(bn, available.classifiers)) {

      if (is(bn, "bn"))
        target = bn$learning$args$training
      else
        target = attr(bn, "training")

    }#THEN
    else {

      stop("missing target node for which to compute the prediction error.")

    }#ELSE

  }#ELSE

  return(target)

}#CHECK.LOSS.TARGET

