\name{graph integration}
\alias{graph integration}
\alias{as.bn.graphNEL}
\alias{as.bn.graphAM}
\alias{as.graphNEL}
\alias{as.graphAM}
\alias{as.graphNEL.bn}
\alias{as.graphNEL.bn.fit}
\alias{as.graphAM.bn}
\alias{as.graphAM.bn.fit}
\title{Import and export networks from the graph package}
\description{

  Convert \code{bn} and \code{bn.fit} objects to \code{graphNEL} and
  \code{graphAM} objects and vice versa.

}
\usage{
\method{as.bn}{graphNEL}(x, ..., check.cycles = TRUE)
\method{as.bn}{graphAM}(x, ..., check.cycles = TRUE)
\method{as.graphNEL}{bn}(x)
\method{as.graphNEL}{bn.fit}(x)
\method{as.graphAM}{bn}(x)
\method{as.graphAM}{bn.fit}(x)
}
\arguments{
  \item{x}{an object of class \code{bn}, \code{bn.fit}, \code{graphNEL},
    \code{graphAM}.}
  \item{\dots}{extra arguments from the generic method (currently ignored).}
  \item{check.cycles}{a boolean value. If \code{FALSE} the returned network will
    not be checked for cycles.}
}
\value{

  An object of the relevant class.

}
\examples{
\dontrun{
library(graph)
a = bn.fit(hc(learning.test), learning.test)
b = as.graphNEL(a)
c = as.bn(b)}}
\author{Marco Scutari}
\keyword{utilities}
\keyword{interfaces to other packages}
