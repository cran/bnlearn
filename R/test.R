
# conditional and unconditional independence tests (vectorized in x, scalar
# in y).
indep.test = function(x, y, sx, data, test, extra.args = list(), alpha = 1,
    learning = TRUE) {

  .Call(call_indep_test,
        x = x,
        y = y,
        sx = sx[sx != ""],
        data = data,
        test = test,
        alpha = alpha,
        extra.args = extra.args,
        learning = learning,
        complete = attr(data, "metadata")$complete.nodes)

}#INDEP.TEST

# test against all possible subsets of the conditioning set (scalar in both x
# and y).
allsubs.test = function(x, y, sx, fixed = character(0), data, test,
    extra.args = list(), alpha = 1, min = 0, max = length(sx), debug = FALSE) {

  .Call(call_allsubs_test,
        x = x,
        y = y,
        sx = c(fixed, sx),
        fixed = fixed,
        data = data,
        test = test,
        alpha = alpha,
        extra.args = extra.args,
        min = as.integer(min),
        max = as.integer(min(max, length(sx))),
        complete = attr(data, "metadata")$complete.nodes,
        debug = debug)

}#ALLSUBS.TEST

# test each variable in turn given the rest as a conditioning set.
roundrobin.test = function(x, z, fixed, data, test, extra.args = list(),
    alpha = 1, debug = FALSE) {

  if (length(z) == 0)
    return(structure(numeric(0), names = character(0)))

  .Call(call_roundrobin_test,
        x = x,
        z = z,
        fixed = fixed,
        data = data,
        test = test,
        alpha = alpha,
        extra.args = extra.args,
        complete = attr(data, "metadata")$complete.nodes,
        debug = debug)

}#ROUNDROBIN.TEST

