
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

