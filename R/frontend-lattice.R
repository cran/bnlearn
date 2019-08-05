
# dispatch qqmath() from lattice for fitted bayesian networks.
bn.fit.qqplot = function(fitted, xlab = "Theoretical Quantiles",
    ylab = "Sample Quantiles", main, ...) {

  if (!is(fitted, c("bn.fit", "bn.fit.gnode", "bn.fit.cgnode")))
    stop("fitted must be an object of class 'bn.fit', 'bn.fit.gnode' and 'bn.fit.cgnode'.")

  if (missing(main)) {

    if (is(fitted, c("bn.fit.gnode", "bn.fit.cgnode")))
      main = paste("Normal Q-Q Plot for Node", fitted$node)
    else
      main = "Normal Q-Q Plot"

  }#THEN

  lattice.gaussian.backend(fitted = fitted, type = "qqplot",
    xlab = xlab, ylab = ylab, main = main, ...)

  invisible(NULL)

}#BN.FIT.QQPLOT

# dispatch histogram() from lattice for fitted bayesian networks.
bn.fit.histogram = function(fitted, density = TRUE, xlab = "Residuals",
    ylab = ifelse(density, "Density", ""), main, ...) {

  if (!is(fitted, c("bn.fit", "bn.fit.gnode", "bn.fit.cgnode")))
    stop("fitted must be an object of class 'bn.fit', 'bn.fit.gnode' and 'bn.fit.cgnode'.")

  if (missing(main)) {

    if (is(fitted, c("bn.fit.gnode", "bn.fit.cgnode")))
      main = paste("Histogram of the Residuals for Node", fitted$node)
    else
      main = "Histogram of the Residuals"

  }#THEN

  lattice.gaussian.backend(fitted = fitted,
    type = ifelse(density, "hist-dens", "hist"),
    xlab = xlab, ylab = ylab, main = main, ...)

  invisible(NULL)

}#BN.FIT.HISTOGRAM

# dispatch xyplot() from lattice for fitted bayesian networks.
bn.fit.xyplot = function(fitted, xlab = "Fitted values", ylab = "Residuals",
    main, ...) {

  if (!is(fitted, c("bn.fit", "bn.fit.gnode", "bn.fit.cgnode")))
    stop("fitted must be an object of class 'bn.fit', 'bn.fit.gnode' and 'bn.fit.cgnode'.")

  if (missing(main)) {

    if (is(fitted, c("bn.fit.gnode", "bn.fit.cgnode")))
      main = paste("Residuals vs Fitted for Node", fitted$node)
    else
      main = "Residuals vs Fitted"

  }#THEN

  lattice.gaussian.backend(fitted = fitted, type = "fitted",
    xlab = xlab, ylab = ylab, main = main, ...)

  invisible(NULL)

}#BN.FIT.XYPLOT

# dispatch barchart() from lattice for fitted bayesian networks.
bn.fit.barchart = function(fitted, xlab = "Probabilities", ylab = "Levels",
    main, ...) {

  if (is(fitted, "bn.fit"))
    stop("only plots of single, discrete nodes are implemented.")
  if (!is(fitted, c("bn.fit.dnode", "bn.fit.onode")))
    stop("fitted must be an object of class 'bn.fit.dnode' or 'bn.fit.onode'.")

  if (missing(main))
    main = paste("Conditional Probabilities for Node", fitted$node)

  lattice.discrete.backend(fitted = fitted, type = "bar",
    xlab = xlab, ylab = ylab, main = main, ...)

  invisible(NULL)

}#BN.FIT.BARCHART

# dispatch dotplot() from lattice for fitted bayesian networks.
bn.fit.dotplot = function(fitted, xlab = "Probabilities", ylab = "Levels",
    main, ...) {

  if (is(fitted, "bn.fit"))
    stop("only plots of single, discrete nodes are implemented.")
  if (!is(fitted, c("bn.fit.dnode", "bn.fit.onode")))
    stop("fitted must be an object of class 'bn.fit.dnode' or 'bn.fit.onode'.")

  if (missing(main))
    main = paste("Conditional Probabilities for Node", fitted$node)

  lattice.discrete.backend(fitted = fitted, type = "dot",
    xlab = xlab, ylab = ylab, main = main, ...)

  invisible(NULL)

}#BN.FIT.DOTPLOT
