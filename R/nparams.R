
# number of parameters for a DAG and an associated data set.
nparams.backend = function(x, data, debug = FALSE) {

  .Call(call_nparams_cgnet,
        graph = x,
        data = .data.frame.column(data, names(x$nodes)),
        debug = debug)

}#NPARAMS.BACKEND

# number of parameters for a fitted network.
nparams.fitted = function(x, debug = FALSE) {

  .Call(call_nparams_fitted,
        bn = x,
        debug = debug)

}#NPARAMS.FITTED

