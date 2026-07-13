
# check the method used to fit the parameters of the network.
check.fitting.method = function(method, data, allowed = available.fits) {

  # check which type of data we are dealing with.
  type = data.type(data)

  if (!missing(method) && !is.null(method)) {

    # check the fitting method.
    check.label(method, choices = allowed, labels = fits.labels,
      argname = "parameter estimator", see = "bn.fit")
    # check that the method is applicable to the data.
    if ((method %in% available.dbn.fits) && (type %!in% discrete.data.types))
      stop("parameter estimator '", method, "' may only be used with discrete data.")
    if ((method %in% c(available.gbn.fits, available.zibn.fits)) &&
        (type != "continuous"))
      stop("parameter estimator '", method, "' may only be used with continuous data.")
    if ((method %in% available.cgbn.fits) && (type != "mixed-cg"))
      stop("parameter estimator '", method, "' may only be used with a mixture of continuous and discrete data.")

    # count data must contain only round numbers.
    if (method %in% available.zibn.fits) {

      all.nonnegative = sapply(data, function(x) all(x >= 0, na.rm = TRUE))
      if (any(!all.nonnegative))
        stop("parameter estimator '", method, "' cannot be used with negative data.")
      all.counts = sapply(data, function(x) all(x == floor(x), na.rm = TRUE))
      if (any(!all.counts))
        warning("parameter estimator '", method, "' is used with non-round count data.")

    }#THEN

    return(method)

  }#THEN
  else {

    if (type %in% discrete.data.types)
      return("mle")
    else if (type == "continuous")
      return("mle-g")
    else if (type == "mixed-cg")
      return("mle-cg")

  }#ELSE

}#CHECK.FITTING.METHOD

# sanitize the extra arguments passed to fitting functions.
check.fitting.args = function(method, network, data, extra.args) {

  if (has.argument(method, "replace.unidentifiable", fits.extra.args)) {

    if (is.null(extra.args[["replace.unidentifiable"]]))
      extra.args[["replace.unidentifiable"]] = FALSE
    else
      check.logical(extra.args[["replace.unidentifiable"]])

  }#THEN

  # check the imaginary sample size.
  if (has.argument(method, "iss", fits.extra.args))
    extra.args[["iss"]] = check.iss(iss = extra.args[["iss"]], network = network)

  # check the imaginary sample size in the dirichlet hierarchical parameter
  # estimators.
  if (has.argument(method, "alpha0", fits.extra.args))
    extra.args[["alpha0"]] =
      check.iss(iss = extra.args[["alpha0"]], network = network)

  # check grouping variable in hierarchical parameter estimators.
  if (has.argument(method, "group", fits.extra.args)){

    # it is sufficient to check that the grouping variable is in the network
    # because the consistency between the network and the data is checked
    # elsewhere.
    check.nodes(extra.args[["group"]], graph = network, max.nodes = 1)

    # by construction, the grouping variable is the only root node in the
    # network: it has arcs pointing to all other nodes.
    if (length(network$nodes[[extra.args[["group"]]]]$parents) != 0)
      stop("the grouping variable should be a root node in the network.")
    if (length(network$nodes[[extra.args[["group"]]]]$children) !=
        (length(network$nodes) - 1))
      stop("not all nodes are children of the grouping variable.")

  }#THEN

  # check the auxiliary parameter of fitting methods in EM-like approaches.
  if (has.argument(method, "fit", fits.extra.args))
    extra.args[["fit"]] =
      check.fitting.method(method = extra.args[["fit"]], data = data,
        allowed = complete.data.fits)
  if (has.argument(method, "fit.args", fits.extra.args))
    extra.args[["fit.args"]] =
      check.fitting.args(method = extra.args[["fit"]],
        network = network, data = data, extra.args[["fit.args"]])

  # check the auxiliary parameter of imputation methods in EM-like approaches.
  if (has.argument(method, "impute", fits.extra.args))
    extra.args[["impute"]] =
      check.imputation.method(method = extra.args[["impute"]], data = data)
  if (has.argument(method, "impute.args", fits.extra.args))
    extra.args[["impute.args"]] =
      check.imputation.extra.args(method = extra.args[["impute"]],
        extra.args = extra.args[["impute.args"]])

  # check the nunber of iterations of EM-like approaches.
  if (has.argument(method, "max.iter", fits.extra.args))
    extra.args[["max.iter"]] =
      check.max.iter(extra.args[["max.iter"]], default = 5)

  # check the threshold for the log-likelihood difference in EM-like approaches.
  if (has.argument(method, "loglik.threshold", fits.extra.args))
    extra.args[["loglik.threshold"]] =
      check.loglik.threshold(extra.args[["loglik.threshold"]])

  # check the threshold for the parameter differences in EM-like approaches.
  if (has.argument(method, "params.threshold", fits.extra.args))
    extra.args[["params.threshold"]] =
      check.params.threshold(extra.args[["params.threshold"]])

  # check the test data for testing convergence in EM-like approaches.
  if (has.argument(method, "newdata", fits.extra.args))
    extra.args[["newdata"]] =
      check.newdata(newdata = extra.args[["newdata"]], network = network,
        data = data, required = FALSE, allow.missing = TRUE)

  # check the network used for initializing EM-like approaches.
  if (has.argument(method, "start", fits.extra.args))
    extra.args[["start"]] =
      check.start.fitted(start = extra.args[["start"]], network = network,
        data = data)

  # make sure that either the thresholds or the maximum number of iterations are
  # set to stop EM-like algorithms.
  if (has.argument(method, "loglik.threshold", fits.extra.args) &&
      has.argument(method, "params.threshold", fits.extra.args) &&
      has.argument(method, "max.iter", fits.extra.args)) {

    if ((extra.args[["loglik.threshold"]] == 0) &&
        (extra.args[["params.threshold"]] == 0) &&
        is.infinite(extra.args[["max.iter"]]))
      stop("'max.iter' cannot be infinite while 'loglik.threshold' and 'params.threshold' are both zero.")

  }#THEN

  # check the controls of the EM estimator for zero-inflated count nodes.
  extra.args = check.em.args(extra.args, fits.extra.args[[method]])

  # warn about and remove unused arguments.
  extra.args = check.unused.args(extra.args, fits.extra.args[[method]])

  return(extra.args)

}#CHECK.FITTING.ARGS

# check the controls of the EM estimator for zero-inflated count nodes: the
# maximum number of EM iterations, the relative convergence tolerance, and the
# M-step type. shared by the bn.fit() methods (mle-zihp, mle-zinb) and by the
# corresponding structure-learning scores; `allowed` is the vector of argument
# names admitted by the method or score, so only relevant controls are checked.
check.em.args = function(extra.args, allowed) {

  # the maximum number of EM iterations and the relative convergence tolerance.
  if ("em.max.iter" %in% allowed)
    extra.args[["em.max.iter"]] =
      check.max.iter(extra.args[["em.max.iter"]], default = 100,
        max = .Machine$integer.max)
  if ("em.tol" %in% allowed) {

    if (is.null(extra.args[["em.tol"]]))
      extra.args[["em.tol"]] = 1e-8
    else if (!is.non.negative(extra.args[["em.tol"]]))
      stop("'em.tol' must be a non-negative numeric value.")

  }#THEN

  # the EM M-step: re-fit each GLM to convergence ("full", the default and the
  # paper-faithful choice) or take a single IRLS step per iteration ("one-step",
  # a cheaper generalized-EM update).
  if ("m.step" %in% allowed) {

    if (is.null(extra.args[["m.step"]]))
      extra.args[["m.step"]] = "full"
    else if (!identical(extra.args[["m.step"]], "full") &&
             !identical(extra.args[["m.step"]], "one-step"))
      stop("'m.step' must be one of \"full\" or \"one-step\".")

  }#THEN

  return(extra.args)

}#CHECK.EM.ARGS

# check the relative log-likelihood threshold for improvement.
check.loglik.threshold = function(threshold) {

  if (!is.null(threshold)) {

    if (!is.non.negative(threshold) && !is.infinite(threshold))
      stop("the log-likelihood threshold must be a non-negative numeric value.")

  }#THEN
  else {

    threshold = 1e-3

  }#ELSE

  return(threshold)

}#CHECK.LOGLIK.THRESHOLD

# check the relative threshold for parameter changes between iterations.
check.params.threshold = function(threshold) {

  if (!is.null(threshold)) {

    if (!is.non.negative(threshold) && !is.infinite(threshold))
      stop("the parameters threshold must be a non-negative numeric value.")

  }#THEN
  else {

    threshold = 1e-3

  }#ELSE

  return(threshold)

}#CHECK.PARAMS.THRESHOLD

# check the fitted network used for initializing EM-like approaches.
check.start.fitted = function(start, network, data) {

  latent = attr(data, "metadata")$latent.nodes

  # having no custom starting model is fine, one will be fitted by EM later.
  if (!is.null(start)) {

    check.fit(start)
    # check the starting network against the data.
    check.fit.vs.data(start, data)

    # if there are any latent variables, they must be connected to at least one
    # non latent-variable for the imputation to work.
    if (any(latent)) {

      for (l in names(which(latent))) {

        nb = c(start[[l]]$parents, start[[l]]$children)

        if ((length(nb) == 0) || (all(nb %in% latent)))
          warning("latent node", l, "is not connected to any observed nodes.")

      }#FOR

    }#THEN

  }#THEN
  else {

    # if there are any latent variables in the data, a starting model is
    # required for the initial imputation.
    if (any(latent))
      stop("a starting model is required to initialize parameter estimation.")

  }#ELSE

  # the starting network is just a device for initializing the EM algorithm, it
  # does not have to have the same arcs as the directed acyclic graph in the
  # parameter learning: no sanitization needed.

  return(start)

}#CHECK.START.FITTED
