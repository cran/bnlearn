
# conditional and unconditional independence tests (vectorized in x, scalar
# in y).
indep.test = function(x, y, sx, data, test, B = 0L, alpha = 1, learning = TRUE,
    complete) {

  .Call(call_indep_test,
        x = x,
        y = y,
        sx = sx[sx != ""],
        data = data,
        test = test,
        B = B,
        alpha = alpha,
        learning = learning,
        complete = complete)

}#INDEP.TEST

# test against all possible subsets of the conditioning set (scalar in both x
# and y).
allsubs.test = function(x, y, sx, fixed = character(0), data, test, B = 0L,
    alpha = 1, min = 0, max = length(sx), complete, debug = FALSE) {

  .Call(call_allsubs_test,
        x = x,
        y = y,
        sx = c(fixed, sx),
        fixed = fixed,
        data = data,
        test = test,
        B = B,
        alpha = alpha,
        min = as.integer(min),
        max = as.integer(min(max, length(sx))),
        complete = complete,
        debug = debug)

}#ALLSUBS.TEST

# test each variable in turn given the rest as a conditioning set.
roundrobin.test = function(x, z, fixed, data, test, B = 0L, alpha = 1,
    complete, debug = FALSE) {

  if (length(z) == 0)
    return(structure(numeric(0), names = character(0)))

  .Call(call_roundrobin_test,
        x = x,
        z = z,
        fixed = fixed,
        data = data,
        test = test,
        B = B,
        alpha = alpha,
        complete = complete,
        debug = debug)

}#ROUNDROBIN.TEST

# Mutual Information (discrete data)
mi.test = function(x, y, ndata, gsquare = TRUE, adjusted = FALSE) {

  .Call(call_mi,
        x = x,
        y = y,
        gsquare = gsquare,
        adjusted = adjusted)

}#MI.TEST

