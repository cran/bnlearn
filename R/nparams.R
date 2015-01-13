
nparams.discrete = function(x, data, real = FALSE, debug = FALSE) {

  # works for both dnode and onode objects, they have the same structure.
  .Call("nparams_dnet",
        graph = x,
        data = minimal.data.frame.column(data, names(x$nodes)),
        real = real,
        debug = debug)

}#NPARAMS.DISCRETE

nparams.gaussian = function(x, debug = FALSE) {

  .Call("nparams_gnet",
        graph = x,
        debug = debug)

}#NPARAMS.GAUSSIAN

nparams.mixedcg = function(x, data, debug = FALSE) {

  .Call("nparams_cgnet",
        graph = x,
        data = minimal.data.frame.column(data, names(x$nodes)),
        debug = debug)

}#NPARAMS.MIXEDCG

nparams.fitted = function(x, debug = FALSE) {

  .Call("fitted_nparams",
        bn = x,
        debug = debug)

}#NPARAMS.FITTED

