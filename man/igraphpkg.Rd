\name{igraph integration}
\alias{igraph integration}
\alias{as.bn.igraph}
\alias{as.igraph}
\alias{as.igraph.bn}
\alias{as.igraph.bn.fit}
\title{Import and export networks from the igraph package}
\description{

  Convert \code{bn} and \code{bn.fit} objects to \code{igraph} and vice versa.

}
\usage{
\method{as.bn}{igraph}(x, ..., check.cycles = TRUE)
\method{as.igraph}{bn}(x, ...)
\method{as.igraph}{bn.fit}(x, ...)
}
\arguments{
  \item{x}{an object of class \code{bn}, \code{bn.fit}, or \code{igraph}.}
  \item{\dots}{extra arguments from the generic method (currently ignored).}
  \item{check.cycles}{a boolean value. If \code{FALSE} the returned network will
    not be checked for cycles.}
}
\value{

  An object of the relevant class.

}
\examples{
\dontrun{
a = bn.fit(hc(learning.test), learning.test)
b = as.igraph(a)
plot(b, edge.arrow.mode = 2L * !igraph::which_mutual(b))
c = as.bn(b)}:}
\author{Marco Scutari}
\keyword{utilities}
\keyword{interfaces to other packages}
