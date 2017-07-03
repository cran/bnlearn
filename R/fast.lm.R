fast.lm = function(data, node, parents, keep.fitted = FALSE) {

  y = minimal.data.frame.column(data, node)

  .Call(call_fast_lm,
        data = data,
        node = node,
        parents = parents,
        keep = keep.fitted)

}#FAST.LM

fast.cglm = function(data, node, parents, configs, keep.fitted = FALSE) {

  .Call(call_fast_cglm,
        data = data,
        node = node,
        parents = parents,
        configs = configs,
        keep = keep.fitted)

}#FAST.CGLM
