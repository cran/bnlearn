
# dispatch qqmath() from lattice for fitted bayesian networks.
bn.fit.qqplot = function(fitted, xlab = "Theoretical Quantiles",
    ylab = "Sample Quantiles", main = "Normal Q-Q Plot", ...) {

  lattice.gaussian.backend(fitted = fitted, type = "qqplot",
    xlab = xlab, ylab = ylab, main = main, ...)

}#BN.FIT.QQPLOT

# dispatch histogram() from lattice for fitted bayesian networks.
bn.fit.histogram = function(fitted, density = TRUE, xlab = "Residuals",
    ylab = ifelse(density, "Density", ""), main = "Histogram of the residuals",
     ...) {

  lattice.gaussian.backend(fitted = fitted,
    type = ifelse(density, "hist-dens", "hist"),
    xlab = xlab, ylab = ylab, main = main, ...)

}#BN.FIT.HISTOGRAM

# dispatch xyplot() from lattice for fitted bayesian networks.
bn.fit.xyplot = function(fitted, xlab = "Fitted values",
   ylab = "Residuals", main = "Residuals vs Fitted", ...) {

  lattice.gaussian.backend(fitted = fitted, type = "fitted",
    xlab = xlab, ylab = ylab, main = main, ...)

}#BN.FIT.XYPLOT

# dispatch barchart() from lattice for fitted bayesian networks.
bn.fit.barchart = function(fitted, xlab = "Probabilities",
   ylab = "Levels", main = "Conditional Probabilities", ...) {

  lattice.discrete.backend(fitted = fitted, type = "bar",
    xlab = xlab, ylab = ylab, main = main, ...)

}#BN.FIT.BARCHART

# dispatch dotplot() from lattice for fitted bayesian networks.
bn.fit.dotplot = function(fitted, xlab = "Probabilities",
   ylab = "Levels", main = "Conditional Probabilities", ...) {

  lattice.discrete.backend(fitted = fitted, type = "dot",
    xlab = xlab, ylab = ylab, main = main, ...)

}#BN.FIT.DOTPLOT
