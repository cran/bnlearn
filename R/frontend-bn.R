
# the generic as method for class bn.
as.bn = function(x, debug = FALSE) {

  UseMethod("as.bn")

}#AS.BN

# get the number of paraters of the bayesian network.
nparams = function(x, data, debug = FALSE) {

  # check x's class.
  check.bn(x)
  # check the data are there.
  check.data(data)
  # check the network against the data
  check.bn.vs.data(x, data)
  # check debug.
  check.logical(debug)
  # nparams is unknown for partially directed graphs.
  if (is.pdag(x$arcs, names(x$nodes)))
    stop("the graph is only partially directed.")

  if (is.data.discrete(data))
    params = nparams.discrete(x, data, real = TRUE)
  else
    params = nparams.gaussian(x)

  if (debug) {

    cat(">  number of parameters per node:\n")
    print(params)

  }#THEN

  return(sum(params))

}#NPARAMS

