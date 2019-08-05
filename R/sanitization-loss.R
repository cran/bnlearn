
# check loss functions' labels.
check.loss = function(loss, data, bn) {

  # check which type of data we are dealing with.
  type = data.type(data)

  if (!is.null(loss)) {

    # check the loss function.
    check.label(loss, choices = loss.functions, labels = loss.labels,
      argname = "loss function", see = "bn.cv")

    if (loss %in% classifiers.loss.functions) {

      if ((is.character(bn) && (bn %!in% classifiers)) ||
          (!is.character(bn) && !is(bn, c("bn.naive", "bn.tan"))))
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

    if ((is.character(bn) && (bn %in% classifiers)) ||
         is(bn, c("bn.naive", "bn.tan")))
      return("pred-exact")
    if (type %in% discrete.data.types)
      return("logl")
    else if (type == "continuous")
      return("logl-g")
    else if (type == "mixed-cg")
      return("logl-cg")

  }#ELSE

}#CHECK.LOSS

# sanitize the extra arguments passed to loss functions.
check.loss.args = function(loss, bn, nodes, data, extra.args) {

  valid.args = loss.extra.args[[loss]]

  if (loss %in% c("pred", "pred-exact", "pred-lw", "pred-lw-cg", "cor",
                  "cor-lw", "cor-lw-cg", "mse", "mse-lw", "mse-lw-cg")) {

    if (!is.null(extra.args$target)) {

      if (!is.string(extra.args$target) || (extra.args$target %!in% nodes))
        stop("target node must be a single, valid node label for the network.")

      # in hybrid networks, check the target has the right data type.
      if (loss %in% c("cor-lw-cg", "mse-lw-cg"))
        if (!is(data[, extra.args$target], "numeric"))
          stop("the target node must be a continuous variable.")
      if (loss == "pred-lw-cg")
        if (!is(data[, extra.args$target], "factor"))
          stop("the target node must be a factor.")

    }#THEN
    else {

      # the target node is obvious for classifiers.
      if (is(bn, c("bn.naive", "bn.tan"))) {

        if (is(bn, "bn"))
          extra.args$target = bn$learning$args$training
        else
          extra.args$target = attr(bn, "training")

      }#THEN
      else {

        stop("missing target node for which to compute the prediction error.")

      }#ELSE

    }#ELSE

    # check the prior distribution.
    if ((is.string(bn) && (bn %in% classifiers)) || is(bn, c("bn.naive", "bn.tan"))) {

      extra.args$prior = check.classifier.prior(extra.args$prior, data[, extra.args$target])
      valid.args = c(valid.args, "prior")

    }#THEN

  }#THEN

  if (loss %in% c("pred-lw", "pred-lw-cg", "cor-lw", "cor-lw-cg", "mse-lw",
                  "mse-lw-cg")) {

    # number of particles for likelihood weighting.
    if (!is.null(extra.args$n)) {

      if (!is.positive.integer(extra.args$n))
        stop("the number of observations to be sampled must be a positive integer number.")

    }#THEN
    else {

      extra.args$n = 500

    }#ELSE

    # which nodes to predict from.
    if (!is.null(extra.args$from))
      check.nodes(extra.args$from, graph = names(data), min.nodes = 1)
    else
      extra.args$from = setdiff(names(data), extra.args$target)

  }#THEN

  # warn about unused arguments.
  check.unused.args(extra.args, valid.args)

  return(extra.args)

}#CHECK.LOSS.ARGS

