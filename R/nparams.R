
# number of parameters for a DAG and an associated data set.
nparams.backend = function(x, data, debug = FALSE) {

  .Call("nparams_cgnet",
        graph = x,
        data = minimal.data.frame.column(data, names(x$nodes)),
        debug = debug)

}#NPARAMS.BACKEND

# number of parameters for a fitted network.
nparams.fitted = function(x, effective = FALSE, debug = FALSE) {

  .Call("nparams_fitted",
        bn = x,
        effective = effective,
        debug = debug)

}#NPARAMS.FITTED

