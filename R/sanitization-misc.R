
# check labels for various arguments.
check.label = function(arg, choices, labels, argname, see) {

  if (!is.string(arg))
    stop("the ", argname, " must be a single character string.")

  if (arg %in% choices)
    return(invisible(NULL))

  # concatenate valid values, optinally with labels.
  if (missing(labels)) {

    choices = paste(paste('"', choices, '"', sep = ""), collapse = ", ")

  }#THEN
  else {

    labels = paste("(", labels[choices], ")", sep = "")
    choices = paste('"', choices, '"', sep = "")
    nl = length(labels)
    choices = paste(choices, labels, collapse = ", ")

  }#THEN

  # mention the most relevant manual page.
  if (missing(see))
    see = character(0)
  else
    see = paste(" See ?", see, " for details.", sep = "")

  # build the error message.
  errmsg = paste("valid ", argname, "(s) are ", choices, ".", see, sep = "")

  # print make sure that it is not truncated if possible at all.
  errlen = unlist(options("warning.length"), use.names = FALSE)
  options("warning.length" = max(1000, min(8170, nchar(errmsg) + 20)))

  stop(errmsg)

  options("warning.length" = errlen)

}#CHECK.LABEL

# check a string that may be a score or a test label.
check.criterion = function(criterion, data) {

  if (!missing(criterion) && !is.null(criterion)) {

    # check and return errors from minimal.check.labels().
    check.label(criterion, choices = c(available.tests, available.scores),
      labels = c(test.labels, score.labels), argname = "criterion",
      see = "bnlearn-package")

  }#THEN
  else {

    # set the defaults using check.score() and check.test().
    if (criterion %in% available.tests)
      criterion = check.test(criterion, data)
    else if (criterion %in% available.scores)
      criterion = check.score(criterion, data)

  }#ELSE

  return(criterion)

}#CHECK.CRITERION

# check the method used for prediction.
check.prediction.method = function(method, data) {

  if (!missing(method) && !is.null(method)) {

    check.label(method, choices = available.prediction.methods,
      labels = prediction.labels, argname = "prediction method", see = "predict")

    return(method)

  }#THEN
  else {

    return("parents")

  }#ELSE

}#CHECK.PREDICTION.METHOD

# check the method used for imputation.
check.imputation.method = function(method, data) {

  if (!missing(method) && !is.null(method)) {

    check.label(method, choices = available.imputation.methods,
      labels = imputation.labels, argname = "imputation method", see = "impute")

    return(method)

  }#THEN
  else {

    return("bayes-lw")

  }#ELSE

}#CHECK.IMPUTATION.METHOD

# check the estimator for the mutual information.
check.mi.estimator = function(estimator, data) {

  # check which type of data we are dealing with.
  type = data.type(data)

  if (!missing(estimator) && !is.null(estimator)) {

    check.label(estimator, choices = available.mi,
      labels = mi.estimator.labels, argname = "mutual information estimator")

    # check if it's the right estimator for the data (discrete, continuous).
    if ((type %!in% discrete.data.types) &&
        (estimator %in% available.discrete.mi))
      stop("estimator '", estimator, "' may be used with discrete data only.")
    if ((type != "continuous") && (estimator %in% available.continuous.mi))
      stop("estimator '", estimator, "' may be used with continuous data only.")

    return(estimator)

  }#THEN
  else {

    if (type %in% discrete.data.types)
      return("mi")
    else
      return("mi-g")

  }#ELSE

}#CHECK.MI.ESTIMATOR

# sanitize the extra arguments passed to the imputation methods.
check.imputation.extra.args = function(method, extra.args) {

  if (method == "parents") {

    # nothing to do.

  }#THEN
  else if (method == "bayes-lw") {

    if (is.null(extra.args$n))
      extra.args$n = 500
    else if (!is.positive.integer(extra.args$n))
      stop("the number of observations to be sampled must be a positive integer number.")

  }#THEN

  check.unused.args(extra.args, imputation.extra.args[[method]])

  return(extra.args)

}#CHECK.IMPUTATION.EXTRA.ARGS

# sanitize the extra arguments passed to the random graph generation algorithms.
check.graph.generation.args = function(method, nodes, extra.args) {

  if (method == "ordered") {

    if (!is.null(extra.args$prob)) {

      # prob must be numeric.
      if (!is.probability(extra.args$prob))
        stop("the branching probability must be a numeric value in [0,1].")

    }#THEN
    else {

      # this default produces graphs with about the same number of
      # arcs as there are nodes.
      extra.args$prob = 2 / (length(nodes) - 1)

    }#ELSE

  }#THEN
  else if (method %in% c("ic-dag", "melancon")) {

    if (!is.null(extra.args$every)) {

      if (!is.positive(extra.args$every))
        stop("'every' must be a positive integer number.")

    }#THEN
    else {

      extra.args$every = 1

    }#ELSE

    if (!is.null(extra.args$burn.in)) {

      if (!is.positive(extra.args$burn.in))
        stop("the burn in length must be a positive integer number.")

    }#THEN
    else {

      extra.args$burn.in = 6 * length(nodes)^2

    }#ELSE

    if (!is.null(extra.args$max.in.degree)) {

      if (!is.positive.integer(extra.args$max.in.degree))
        stop("the maximum in-degree must be a positive integer number.")

      if (extra.args$max.in.degree >= length(nodes)) {

        warning("a node cannot have an in-degree greater or equal to the number of nodes in the graph.")
        warning("the condition on the in-degree will be ignored.")

      }#THEN

    }#THEN
    else {

      extra.args$max.in.degree = Inf

    }#ELSE

    if (!is.null(extra.args$max.out.degree)) {

      if (!is.positive.integer(extra.args$max.out.degree))
        stop("the maximum out-degree must be a positive integer number.")

      if (extra.args$max.out.degree >= length(nodes)) {

        warning("a node cannot have an out-degree greater or equal to the number of nodes in the graph.")
        warning("the condition on the out-degree will be ignored.")

      }#THEN

    }#THEN
    else {

      extra.args$max.out.degree = Inf

    }#ELSE

    if (!is.null(extra.args$max.degree)) {

      if (!is.positive.integer(extra.args$max.degree))
        stop("the maximum out-degree must be a positive integer number.")

      if (is.finite(extra.args$max.in.degree) &&
          extra.args$max.in.degree > extra.args$max.degree)
        stop("the maximun in-degree must be lesser or equal to the maximum degree.")

      if (is.finite(extra.args$max.out.degree) &&
          extra.args$max.out.degree > extra.args$max.degree)
        stop("the maximun out-degree must be lesser or equal to the maximum degree.")

      if (extra.args$max.degree >= length(nodes)) {

        warning("a node cannot have a degree greater or equal to the number of nodes in the graph.")
        warning("the condition on the degree will be ignored.")

      }#THEN

    }#THEN
    else {

      extra.args$max.degree = Inf

    }#ELSE

  }#THEN

  check.unused.args(extra.args, graph.generation.extra.args[[method]])

  return(extra.args)

}#CHECK.GRAPH.GENERATION.ARGS

# sanitize the extra arguments passed to Bayesian classifiers.
check.classifier.args = function(method, data, training, explanatory,
    extra.args) {

  if (method == "tree.bayes") {

    # check the label of the mutual information estimator.
    extra.args$estimator = check.mi.estimator(extra.args$estimator, data)

    # check the node to use the root of the tree (if not specified pick the first
    # explanatory variable assuming natural ordering).
    if (!is.null(extra.args$root))
      check.nodes(extra.args$root, graph = explanatory, max.nodes = 1)
    else
      extra.args$root = explanatory[1]

  }#THEN

  return(extra.args)

}#CHECK.CLASSIFIER.ARGS

# warn about unused arguments.
check.unused.args = function(dots, used.args) {

  if (is(dots, "list"))
    args = names(dots)
  else
    args = dots

  unused = args %!in% used.args

  if (any(unused))
    warning("unused argument(s):", paste0(" '", args[unused], "'"), ".")

}#CHECK.UNUSED.ARGS

check.amat = function(amat, nodes) {

  # a node is needed.
  if (missing(amat))
    stop("no adjacency matrix specified.")
  # the adjacency matrix must, well, be a matrix.
  if (!is(amat, "matrix") || (ncol(amat) != nrow(amat)) || (length(dim(amat)) != 2))
    stop("an adjacency matrix must be a 2-dimensional square matrix.")
  # check the dimensions against the number of nodes in the graph.
  if (any(dim(amat) != length(nodes)))
    stop("the dimensions of the adjacency matrix do not agree with the number of nodes in the graph.")
  # column names must be valid node labels.
  if (!is.null(colnames(amat)))
    if (any(colnames(amat) %!in% nodes))
      stop("node (column label) not present in the graph.")
  # column names must be valid node labels.
  if (!is.null(rownames(amat)))
    if (any(rownames(amat) %!in% nodes))
      stop("node (row label) not present in the graph.")
  # column names must match with row names.
  if (!is.null(colnames(amat)) && !is.null(rownames(amat))) {

    if (!identical(colnames(amat), rownames(amat)))
      stop("row/column names mismatch in the adjacency matrix.")

    if (!identical(colnames(amat), nodes) || !identical(rownames(amat), nodes)) {

      warning("rearranging the rows/columns of the adjacency matrix.")

      amat = amat[nodes, nodes, drop = FALSE]

    }#THEN

  }#THEN
  # make really sure the adjacency matrix is made up of integers.
  if (storage.mode(amat) != "integer")
    storage.mode(amat) = "integer"
  # check the elements of the matrix.
  if (!all((amat == 0L) | (amat == 1L)))
    stop("all the elements of an adjacency matrix must be equal to either 0 or 1.")
  # no arcs from a node to itself.
  if (any(diag(amat) != 0))
    stop("the elements on the diagonal must be zero.")

  return(amat)

}#CHECK.AMAT

# check a prior distribution against the observed variable.
check.classifier.prior = function(prior, training) {

  if (missing(prior) || is.null(prior)) {

    # use the empirical probabilities in the fitted network, or a flat prior
    # as a last resort.
    if (is(training, c("bn.fit.dnode", "bn.fit.onode")))
      prior = training$prob
    else
      prior = rep(1, nlevels(training))

  }#THEN
  else {

    if (is(training, c("bn.fit.dnode", "bn.fit.onode")))
      nlvls = dim(training$prob)[1]
    else
      nlvls = nlevels(training)

    if (length(prior) != nlvls)
      stop("the prior distribution and the training variable have a different number of levels.")
    if (!is.nonnegative.vector(prior))
      stop("the prior distribution must be expressed as a probability vector.")

    # make sure the prior probabilities sum to one.
    prior = prior / sum(prior)

  }#ELSE

  return(prior)

}#CHECK.CLASSIFIER.PRIOR

# check a vector of weights.
check.weights = function(weights, len) {

  if (missing(weights) || is.null(weights)) {

    weights = rep(1, len)

  }#THEN
  else {

    if (!is.nonnegative.vector(weights))
      stop("missing or negative weights are not allowed.")

    if (length(weights) != len)
      stop("wrong number of weights, ", length(weights),
        " weights while ", len, " are needed.")

    weights = prop.table(weights)

  }#ELSE

  return(weights)

}#CHECK.WEIGHTS

