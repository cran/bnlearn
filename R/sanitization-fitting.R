
# check the method used to fit the parameters of the network.
check.fitting.method = function(method, data) {

  # check which type of data we are dealing with.
  type = data.type(data)

  if (!missing(method) && !is.null(method)) {

    # check the fitting method.
    check.label(method, choices = available.fits, labels = fits.labels,
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

    if (is.null(extra.args$replace.unidentifiable))
      extra.args$replace.unidentifiable = FALSE
    else
      check.logical(extra.args$replace.unidentifiable)

  }#THEN

  # check the imaginary sample size.
  if (has.argument(method, "iss", fits.extra.args))
    extra.args$iss = check.iss(iss = extra.args$iss, network = network)

  # check the imaginary sample size in the dirichlet hierarchical parameter
  # estimators.
  if (has.argument(method, "alpha0", fits.extra.args))
    extra.args$alpha0 = check.iss(iss = extra.args$alpha0, network = network)

  # check grouping variable in hierarchical parameter estimators.
  if (has.argument(method, "group", fits.extra.args)){

    # it is sufficient to check that the grouping variable is in the network
    # because the consistency between the network and the data is checked
    # elsewhere.
    check.nodes(extra.args$group, graph = network, max.nodes = 1)

    # by construction, the grouping variable is the only root node in the
    # network: it has arcs pointing to all other nodes.
    if (length(network$nodes[[extra.args$group]]$parents) != 0)
      stop("the grouping variable should be a root node in the network.")
    if (length(network$nodes[[extra.args$group]]$children) !=
        (length(network$nodes) - 1))
      stop("not all nodes are children of the grouping variable.")

  }#THEN

  # warn about unused arguments.
  check.unused.args(extra.args, fits.extra.args[[method]])

  return(extra.args)

}#CHECK.FITTING.ARGS

