\name{plot.bn}
\alias{plot.bn}
\title{Plot a Bayesian network}
\description{

  Plot the graph associated with a small Bayesian network.

}
\usage{
\method{plot}{bn}(x, ylim = c(0,600), xlim = ylim, radius = 250,
  arrow = 35, highlight = NULL, color = "red", ...)
}
\arguments{
  \item{x}{an object of class \code{bn}.}
  \item{ylim}{a numeric vector with two components containing the range of the
    y-axis.}
  \item{xlim}{a numeric vector with two components containing the range of the
    x-axis.}
  \item{radius}{a numeric value containing the radius of the nodes.}
  \item{arrow}{a numeric value containing the length of the arrow heads.}
  \item{highlight}{a vector of character strings, representing the labels of
    the nodes (and corresponding arcs) to be highlighted.}
  \item{color}{an integer or character string (the highlight colour).}
  \item{\dots}{other graphical parameters to be passed through to plotting
    functions.}
}
\note{
  The following arguments are always overridden:
  \itemize{
    \item \code{axes} is set to \code{FALSE}.
    \item \code{xlab} is set to an empty string.
    \item \code{ylab} is set to an empty string.
  }
}
\examples{
data(learning.test)
cpdag = pc.stable(learning.test)

plot(cpdag)

## highlight node B and related arcs.
plot(cpdag, highlight = "B")
## highlight B and its Markov blanket.
plot(cpdag, highlight = c("B", mb(cpdag, "B")))

## a more compact plot.
par(oma = rep(0, 4), mar = rep(0, 4), mai = rep(0, 4),
  plt = c(0.06, 0.94, 0.12, 0.88))
plot(cpdag)
}
\seealso{\code{\link{graphviz.plot}}.}
\author{Marco Scutari}
\keyword{plots}
