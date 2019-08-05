
# check the method used to fit the parameters of the network.
check.fitting.method = function(method, data) {

  # check which type of data we are dealing with.
  type = data.type(data)

  if (!is.null(method)) {

    # check the fitting method.
    check.label(method, choices = available.fitting.methods,
      labels = fitting.labels, argname = "fitting method", see = "bn.fit")
    # Bayesian parameter estimation is implemented only for discrete data.
    if ((type %in% c("continuous", "mixed-cg")) && (method == "bayes"))
      stop("Bayesian parameter estimation for (conditional) Gaussian Bayesian networks is not implemented.")

    return(method)

  }#THEN
  else {

    return("mle")

  }#ELSE

}#CHECK.FITTING.METHOD

# sanitize the extra arguments passed to fitting functions.
check.fitting.args = function(method, network, data, extra.args) {

  if (has.argument(method, "replace.unidentifiable", fitting.extra.args)) {

    if (is.null(extra.args$replace.unidentifiable))
      extra.args$replace.unidentifiable = FALSE
    else
      check.logical(extra.args$replace.unidentifiable)

  }#THEN

  # check the imaginary sample size.
  if (has.argument(method, "iss", fitting.extra.args))
    extra.args$iss = check.iss(iss = extra.args$iss, network = network)

  # warn about unused arguments.
  check.unused.args(extra.args, fitting.extra.args[[method]])

  return(extra.args)

}#CHECK.FITTING.ARGS

