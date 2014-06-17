
# conditional and unconditional independence tests (vectorized in x).
indep.test = function(x, y, sx, data, test, B = 0, alpha = 1, learning = TRUE) {

  .Call("indep_test",
        x = x,
        y = y,
        sx = sx[sx != ""],
        data = data,
        test = test,
        B = B,
        alpha = alpha,
        learning = learning)

}#INDEP.TEST

# test against all possible subsets of the conditioning set.
allsubs.test = function(x, y, sx, data, test, B = 0, alpha = 1, min = 0,
    max = length(sx), debug = FALSE) {

  sx = sx[sx != ""]
  if (missing(max))
    max = length(sx)

  .Call("allsubs_test",
        x = x,
        y = y,
        sx = sx,
        data = data,
        test = test,
        B = B,
        alpha = alpha,
        min = as.integer(min),
        max = as.integer(max),
        debug = debug)

}#ALLSUBS.TEST

# Mutual Information (discrete data)
mi.test = function(x, y, ndata, gsquare = TRUE, adjusted = FALSE) {

  .Call("mi",
        x = x,
        y = y,
        gsquare = gsquare,
        adjusted = adjusted)

}#MI.TEST

