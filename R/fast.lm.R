# fast implementation of linear models.
fast.lm = function(data, node, parents, keep.fitted = FALSE,
     with.missing = FALSE) {

  .Call(call_fast_lm,
        data = data,
        node = node,
        parents = parents,
        keep = keep.fitted,
        missing = with.missing)

}#FAST.LM

# fast implementation of mixtures of linear models.
fast.cglm = function(data, node, parents, configs, keep.fitted = FALSE,
    with.missing = FALSE) {

  .Call(call_fast_cglm,
        data = data,
        node = node,
        parents = parents,
        configs = configs,
        keep = keep.fitted,
        missing = with.missing)

}#FAST.CGLM
