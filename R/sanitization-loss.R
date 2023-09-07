
# check loss functions' labels.
check.loss = function(loss, data, bn) {

  # check which type of data we are dealing with.
  type = data.type(data)

  if (!is.null(loss)) {

    # check the loss function.
    check.label(loss, choices = loss.functions, labels = loss.labels,
      argname = "loss function", see = "bn.cv")

    if (loss %in% classifiers.loss.functions) {

      if ((is.character(bn) && (bn %!in% classification.algorithms)) ||
          (!is.character(bn) && !is(bn, available.classifiers)))
        stop("loss function '", loss, "' may only be used with classifiers.")

    }#THEN

    if ((type %!in% discrete.data.types) && (loss %in% discrete.loss.functions))
      stop("loss function '", loss, "' may only be used with discrete data.")
    if ((type != "continuous") && (loss %in% continuous.loss.functions))
      stop("loss function '", loss, "' may only be used with continuous data.")
    if ((type != "mixed-cg") && (loss %in% mixedcg.loss.functions))
      stop("loss function '", loss, "' may only be used with a mixture of continuous and discrete data.")

    return(loss)

  }#THEN
  else {

    if ((is.character(bn) && (bn %in% classification.algorithms)) ||
         is(bn, available.classifiers))
      return("pred-exact")
    if (type %in% discrete.data.types)
      return("logl")
    else if (type == "continuous")
      return("logl-g")
    else if (type == "mixed-cg")
      return("logl-cg")

  }#ELSE

}#CHECK.LOSS

# check the target node the loss is computed for.
check.loss.target = function(target, loss, data, bn, nodes) {

  if (!is.null(target)) {

    if (!is.string(target) || (target %!in% nodes))
      stop("target node must be a single, valid node label for the network.")

    # in hybrid networks, check the target has the right data type.
    if (loss %in% c("cor-lw-cg", "mse-lw-cg"))
      if (!is(data[, target], "numeric"))
        stop("the target node must be a continuous variable.")
    if (loss == "pred-lw-cg")
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

# sanitize the extra arguments passed to loss functions.
check.loss.args = function(loss, bn, nodes, data, extra.args) {

  # target node that that the loss is computed for.
  if (has.argument(loss, "target", loss.extra.args))
    extra.args[["target"]] =
      check.loss.target(target = extra.args[["target"]], loss = loss,
        data = data, bn = bn, nodes = nodes)

  # check the prior distribution for the target variable of a classifier.
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

