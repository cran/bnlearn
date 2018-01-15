
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

