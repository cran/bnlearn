
dedup.backend = function(data, threshold, complete, debug = FALSE) {

  .Call(call_dedup,
        data = data,
        threshold = threshold,
        complete = complete,
        debug = debug)

}#DEDUP.BACKEND

discretize.backend = function(data, method, breaks, ordered = FALSE, extra.args,
    debug = FALSE) {

  if (method %in% c("quantile", "interval")) {

    marginal.discretize.backend(data = data, method = method, breaks = breaks,
      ordered = ordered, debug = debug)

  }#THEN
  else if (method == "hartemink") {

    joint.discretize.backend(data = data, method = method, breaks = breaks,
      ordered = ordered, initial.discretization = extra.args$idisc,
      initial.breaks = extra.args$ibreaks, debug = debug)

  }#ELSE

}#DISCRETIZE.BACKEND

marginal.discretize.backend = function(data, method, breaks, ordered = FALSE,
    debug = FALSE) {

  .Call(call_marginal_discretize,
        data = data,
        method = method,
        breaks = as.integer(breaks),
        ordered = ordered,
        debug = debug)

}#MARGINAL.DISCRETIZE.BACKEND

joint.discretize.backend = function(data, method, breaks, ordered = FALSE,
    initial.discretization = "quantile", initial.breaks, debug = debug) {

  .Call(call_joint_discretize,
        data = data,
        method = method,
        breaks = as.integer(breaks),
        ordered = ordered,
        initial.discretization = initial.discretization,
        initial.breaks = as.integer(initial.breaks),
        debug = debug)

}#JOINT.DISCRETIZE.BACKEND
