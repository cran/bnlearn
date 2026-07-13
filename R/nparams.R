
# number of parameters for a DAG and an associated data set.
nparams.backend = function(dag, data, estimator, debug = FALSE) {

  .Call(call_nparams_structure,
        graph = dag,
        data = data[, names(dag$nodes), drop = TRUE],
        estimator = estimator,
        debug = debug)

}#NPARAMS.BACKEND

# number of parameters for a fitted network.
nparams.fitted = function(bn, debug = FALSE) {

  .Call(call_nparams_fitted,
        bn = bn,
        debug = debug)

}#NPARAMS.FITTED

