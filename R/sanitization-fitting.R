
# check the method used to fit the parameters of the network.
check.fitting.method = function(method, data, allowed = available.fits) {

  # check which type of data we are dealing with.
  type = data.type(data)

  if (!missing(method) && !is.null(method)) {

    # check the fitting method.
    check.label(method, choices = allowed, labels = fits.labels,
      argname = "fitting method", see = "bn.fit")
    # check that the method is applicable to the data.
    if ((method %in% available.dbn.fits) && (type %!in% discrete.data.types))
      stop("parameter estimator '", method, "' may only be used with discrete data.")
    if ((method %in% available.gbn.fits) && (type != "continuous"))
      stop("parameter estimator '", method, "' may only be used with continuous data.")
    if ((method %in% available.cgbn.fits) && (type != "mixed-cg"))
      stop("paramter estimator '", method, "' may only be used with a mixture of continuous and discrete data.")

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
  if (has.argument(method, "threshold", fits.extra.args))
    extra.args[["threshold"]] =
      check.loglik.threshold(extra.args[["threshold"]])

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

  # make sure that the threshold and the maximum number of iterations of EM-like
  # algorithms cannot be infinite at the same time.
  if (has.argument(method, "threshold", fits.extra.args) &&
      has.argument(method, "max.iter", fits.extra.args)) {

    if (is.infinite(extra.args[["threshold"]]) &&
        is.infinite(extra.args[["max.iter"]]))
      stop("'threshold' and 'max.iter' cannot be infinite at the same time.")

  }#THEN

  # warn about and remove unused arguments.
  extra.args = check.unused.args(extra.args, fits.extra.args[[method]])

  return(extra.args)

}#CHECK.FITTING.ARGS

# check the relative log-likelihood threshold for improvement.
check.loglik.threshold = function(threshold) {

  if (!is.null(threshold)) {

    if (!is.positive(threshold) && !isTRUE(all.equal(threshold, Inf)))
      stop("the threshold must be a positive numeric value.")

  }#THEN
  else {

    threshold = 1e-3

  }#ELSE

  return(threshold)

}#CHECK.LOGLIK.THRESHOLD

# check the fitted network used for initializing EM-like approaches.
check.start.fitted = function(start, network, data) {

  # having no custom starting model is fine, one will be fitted by EM later.
  if (!is.null(start)) {

    # check the class.
    check.fit(start)
    # check the starting network against the data.
    check.fit.vs.data(start, data)

  }#THEN
  else {

    # if there are any latent variables in the data, a starting model is
    # required for the initial imputation.
    if (any(attr(data, "metadata")$latent.nodes))
      stop("a starting model is required to initialize parameter estimation.")

  }#ELSE

  # the starting network is just a device for initializing the EM algorithm, it
  # does not have to have the same arcs as the directed acyclic graph in the
  # parameter learning: no sanitization needed.

  return(start)

}#CHECK.START.FITTED
