
# compute the log-likelihood of some data for a given network.
loglikelihood = function(fitted, data, by.sample = FALSE, keep = names(fitted),
    propagate.missing = FALSE, as.loss = FALSE, debug = FALSE) {

  .Call(call_loglikelihood_function,
        fitted = fitted,
        data = data,
        by.sample = by.sample,
        keep.nodes = keep,
        propagate.missing = propagate.missing,
        as.loss = as.loss,
        debug = debug)

}#LOGLIKELIHOOD
