\name{gRain integration}
\alias{gRain integration}
\alias{as.bn.fit}
\alias{as.bn.fit.grain}
\alias{as.grain}
\alias{as.grain.bn.fit}
\title{ Import and export networks from the gRain package }
\description{

  Convert \code{bn.fit} objects to \code{grain} objects and vice versa.

}
\usage{
\method{as.bn.fit}{grain}(x)
\method{as.grain}{bn.fit}(x)
}
\arguments{
  \item{x}{an object of class \code{grain} (for \code{as.bn.fit}) or 
    \code{bn.fit} (for \code{as.grain}).}
}
\note{

  Conditional probability tables in \code{grain} objects must be completely
  specified; on the other hand, \code{bn.fit} allows \code{NaN} values for
  unobserved parents' configurations. Such \code{bn.fit} objects cannot be
  converted to \code{grain} objects as they are. A simple solution to this
  problem is to fit another \code{bn.fit} with \code{method = "bayes"} and
  a low \code{iss} value, using the same data and network structure.

}
\value{

  An object of class \code{grain} (for \code{as.grain}) or \code{bn.fit}
  (for \code{as.bn.fit}).

}
\examples{
\dontrun{
library(gRain)
a = bn.fit(hc(learning.test), learning.test)
b = as.grain(a)
c = as.bn.fit(b)}}
\author{ Marco Scutari }
\keyword{utilities}
\keyword{multivariate}
