
nparams.backend = function(x, data, debug = FALSE) {

  # handles all of discrete, Gaussian and conditional Gaussian networks.
  .Call("nparams_cgnet",
        graph = x,
        data = minimal.data.frame.column(data, names(x$nodes)),
        debug = debug)

}#NPARAMS.BACKEND

nparams.fitted = function(x, effective = FALSE, debug = FALSE) {

  .Call("nparams_fitted",
        bn = x,
        effective = effective,
        debug = debug)

}#NPARAMS.FITTED

