
# check bn metadata against the data it's used with.
check.bn.vs.data = function(bn, data) {

  # check which type of data we are dealing with.
  type = data.type(data)

  # the number of variables must be the same
  if (length(names(bn$nodes)) != ncol(data))
    stop("the network and the data have different numbers of variables.")
  # the variables must be the same.
  if (length(setdiff(names(bn$nodes), names(data))) != 0)
    stop("the variables in the data and in the network do not match.")
  # data type versus network structure.
  if (type == "mixed-cg")
    check.arcs.against.assumptions(bn$arcs, data, "mi-cg")

}#CHECK.BN.VS.DATA

# check bn.fit metadata against the data it's used with.
check.fit.vs.data = function(fitted, data, subset) {

  fitted.names = names(fitted)
  # check which type of data we are dealing with.
  dtype = data.type(data)

  if (missing(subset)) {

    # the number of variables must be the same.
    if (length(fitted.names) != ncol(data))
      stop("the network and the data have different numbers of variables.")
    # the variables must be the same.
    if (length(setdiff(fitted.names , names(data))) != 0)
      stop("the variables in the data and in the network do not match.")

    subset = fitted.names

  }#THEN
  else {

    # the number of variables must not exceed that of the network.
    if (length(subset) > length(fitted.names))
      stop("the data have more variables than the network.")
    # all the variables in the subset must be present in the data.
    absent = (subset %!in% names(data))
    if (any(absent))
      stop("required variables '", paste(subset[absent], collapse = " "),
           "' are not present in the data.")
    # all the variables in the subset must also be present in the network.
    absent = (subset %!in% fitted.names)
    if (any(absent))
      stop("required variables '", paste(subset[absent], collapse = " "),
           "' are not present in the network.")

  }#ELSE

  .Call(call_fitted_vs_data,
        fitted = fitted,
        data = data,
        subset = subset)

}#CHECK.FIT.VS.DATA

# check bn.fit.{d,g}node metadata against the data it's used with.
check.fit.node.vs.data = function(fitted, data) {

  relevant = c(fitted$node, fitted$parents)
  # check which type of data we are dealing with.
  type = data.type(data)

  # check whether all relevant nodes are in the data.
  if (any(relevant %!in% names(data)))
    stop("not all required nodes are present in the data.")
  # data type versus network type.
  if (is(fitted, "bn.fit.dnode") && (type == "continuous"))
      stop("continuous data and discrete network.")
  if (is(fitted, "bn.fit.gnode") &&
      (type %in% discrete.data.types))
    stop("discrete data and continuous network.")
  # double-check the levels of the variables against those of the nodes.
  if (is(fitted, "bn.fit.dnode")) {

    for (node in relevant) {

      data.levels = levels(data[, node])
      if (length(relevant) == 1)
        node.levels = dimnames(fitted$prob)[[1]]
      else
        node.levels = dimnames(fitted$prob)[[node]]

      if (!identical(data.levels, node.levels))
        stop("the levels of node '", node, "' do not match the levels of the ",
             "corresponding variable in the data.")

    }#FOR

  }#THEN

}#CHECK.FIT.NODE.VS.DATA

