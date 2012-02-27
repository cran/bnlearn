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
\title{ Import and export networks from the graph package }
\description{

  Convert \code{bn} and \code{bn.fit} objects to \code{graphNEL} and
  \code{graphAM} objects and vice versa.

}
\usage{
\method{as.bn}{graphNEL}(x)
\method{as.bn}{graphAM}(x)
\method{as.graphNEL}{bn}(x)
\method{as.graphNEL}{bn.fit}(x)
\method{as.graphAM}{bn}(x)
\method{as.graphAM}{bn.fit}(x)
}
\arguments{
  \item{x}{an object of class \code{bn}, \code{bn.fit}, \code{graphNEL},
    \code{graphAM}.}
}
\note{

  The corresponding S4 methods are exported as well, and are just wrappers
  around the S3 ones. So, for example, both \code{as.graphNEL(x)} and
  \code{as(x, "graphNEL")} work and return identical objects.

}
\value{

  An object of the relevant class.

}
\examples{
\dontrun{
library(graph)
a = bn.fit(hc(learning.test), learning.test)
b = as.graphNEL(a)
c = as.bn(b)
identical(as.graphNEL(a), as(a, "graphNEL"))}}
\author{ Marco Scutari }
\keyword{utilities}
\keyword{multivariate}
