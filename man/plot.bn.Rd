\name{plot.bn}
\alias{plot.bn}
\title{ Plot a Bayesian network }
\description{

  Plot the graph associated with a small Bayesian network.

}
\usage{
\method{plot}{bn}(x, ylim = c(0,600), xlim = ylim, radius = 250,
  arrow = 35, highlight = NULL, color = "red", ...)
}
\arguments{
  \item{x}{an object of class \code{bn}.}
  \item{ylim}{a numeric vector with two components containing the range
          on y-axis.}
  \item{xlim}{a numeric vector with two components containing the range
          on x-axis.}
  \item{radius}{a numeric value containing the radius of the nodes.}
  \item{arrow}{a numeric value containing the length of the arrow heads.}
  \item{highlight}{a vector of character strings, representing the
    labels of the nodes (and corresponding arcs) to be highlighted.}
  \item{color}{an integer or character string (the highlight colour).}
  \item{\dots}{other parameters to be passed through to plotting
    functions.}
}
\note{
  The following graphical parameters are always overridden:
  \itemize{
    \item \code{axes} is set to \code{FALSE}.
    \item \code{xlab} is set to an empty string.
    \item \code{ylab} is set to an empty string.
  }
}
\examples{
data(learning.test)
res = gs(learning.test)

plot(res)

## highlight node B and related arcs.
plot(res, highlight = "B")
## highlight B and its Markov blanket.
plot(res, highlight = c("B", mb(res, "B")))

## a more compact plot.
par(oma = rep(0, 4), mar = rep(0, 4), mai = rep(0, 4),
  plt = c(0.06, 0.94, 0.12, 0.88))
plot(res)
}
\seealso{\code{\link{graphviz.plot}}.}
\author{ Marco Scutari }
\keyword{hplot}
