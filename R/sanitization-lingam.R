
# check the method label for DirectLiNGAM.
check.dlingam.maximize = function(method) {

  if (missing(method))
    maximize = "alasso"
  else if (!identical(method, "alasso"))
    stop("the 'maximize' learning approach can only be 'alasso'.")

  return(method)

}#CHECK.DLINGAM.MAXIMIZE

# check the optional arguments to DirectLiNGAM second step.
check.dlingam.maximize.args = function(method, extra.args, data) {

  if (method == "alasso") {

    # the exponent of the weights in the adaptive LASSO.
    if (is.null(extra.args$gamma))
      extra.args$gamma = 1
    else if (!is.positive(extra.args$gamma))
      stop("'gamma' must be a positive number, usually 0.5, 1, or 2.")

    # lower threshold lambda values in proportion to the maximum lambda.
    if (is.null(extra.args$lambda.min.ratio)) {

      if (nrow(data) > ncol(data))
        extra.args$lambda.min.ratio = 0.0001
      else
        extra.args$lambda.min.ratio = 0.01

    }#THEN
    else if (!is.probability(extra.args$lambda.min.ratio) ||
             (extra.args$lambda.min.ratio == 1)) {

      stop("'lambda.min.ratio' must be a number in [0, 1).")

    }#THEN

    # bound the number of parents for a single node.
    if (is.null(extra.args$pmax))
      extra.args$pmax = ncol(data) + 1
    else if (!is.non.negative.integer(extra.args$pmax) ||
             (extra.args$pmax > ncol(data) + 1))
      stop("'pmax' must be a positive integer number, at most", ncol(data) + 1, ".")

    # likelihood penalty, defaults to BIC.
    extra.args$k =
      check.penalty(extra.args$k, network = NULL, data = data, score = "bic-g")

    extra.args =
      check.unused.args(extra.args, c("gamma", "lambda.min.ratio", "pmax", "k"))

  }#THEN

  return(extra.args)

}#CHECK.DLINGAM.MAXIMIZE.ARGS
