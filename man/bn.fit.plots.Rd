\name{bn.fit plots}
\alias{bn.fit plots}
\alias{bn.fit.qqplot}
\alias{bn.fit.xyplot}
\alias{bn.fit.histogram}
\alias{bn.fit.barchart}
\alias{bn.fit.dotplot}

\title{ Plot fitted Bayesian networks }
\description{

  Plot functions for the \code{bn.fit}, \code{bn.fit.dnode} and
  \code{bn.fit.gnode} classes, based on the \pkg{lattice} package.

}
\usage{
## for Gaussian Bayesian networks.
bn.fit.qqplot(fitted, xlab = "Theoretical Quantiles",
  ylab = "Sample Quantiles", main = "Normal Q-Q Plot", ...)
bn.fit.histogram(fitted, density = TRUE, xlab = "Residuals",
  ylab = ifelse(density, "Density", ""),
  main = "Histogram of the residuals", ...)
bn.fit.xyplot(fitted, xlab = "Fitted values",
  ylab = "Residuals", main = "Residuals vs Fitted", ...)
## for discrete Bayesian networks
bn.fit.barchart(fitted, xlab = "Probabilities",
  ylab = "Levels", main = "Conditional Probabilities", ...)
bn.fit.dotplot(fitted, xlab = "Probabilities",
  ylab = "Levels", main = "Conditional Probabilities", ...)
}
\arguments{
  \item{fitted}{an object of class \code{bn.fit}, \code{bn.fit.dnode}
    or \code{bn.fit.gnode}.}
  \item{xlab, ylab, main}{the label of the x axis, of the y axis,
    and the plot title.}
  \item{density}{a boolean value. If \code{TRUE} the histogram is
    plotted using relative frequencies, and the matching normal
    density is added to the plot.}
  \item{\dots}{additional arguments to be passed to \pkg{lattice}
   functions.}
}
\details{

  \code{bn.fit.qqplot} draws a quantile-quantile plot of the
    residuals.

  \code{bn.fit.histogram} draws a histogram of the residuals,
    using either absolute or relative frequencies.

  \code{bn.fit.xyplot} plots the residuals versus the fitted
    values.

  \code{bn.fit.barchart} and \code{bn.fit.dotplot} plot
    the probabilities in the conditional probability table
    associated with each node.

}
\value{

  The \pkg{lattice} plot objects. Note that if auto-printing is
  turned off (for example when the code is loaded with the
  \code{source} function), the return value must be printed
  explicitly for the plot to be displayed.

}
\seealso{
  \code{\link{bn.fit}}, \code{\link{bn.fit class}}.
}
\author{ Marco Scutari }
\keyword{hplot}
