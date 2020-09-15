\name{gRain integration}
\alias{gRain integration}
\alias{as.bn.fit}
\alias{as.bn.fit.grain}
\alias{as.bn.grain}
\alias{as.grain}
\alias{as.grain.bn.fit}
\alias{as.grain.bn}
\title{Import and export networks from the gRain package}
\description{

  Convert \code{bn.fit} objects to \code{grain} objects and vice versa.

}
\usage{
\method{as.bn.fit}{grain}(x, including.evidence = FALSE, ...)
\method{as.grain}{bn.fit}(x)
\method{as.bn}{grain}(x, ...)
}
\arguments{
  \item{x}{an object of class \code{grain(code)} (for \code{as.bn.fit}) or
    \code{bn.fit()} (for \code{as.grain}).}
  \item{including.evidence}{a boolean value. If \code{FALSE}, the \code{grain}
    object is converted without considering any evidence that has been set into
    it. If \code{TRUE}, any hard evidence is carried over into the \code{bn.fit}
    object as a zero-one probability distribution.}
  \item{\dots}{extra arguments from the generic method (currently ignored).}
}
\note{

  Conditional probability tables in \code{grain} objects must be completely
  specified; on the other hand, \code{bn.fit} allows \code{NaN} values for
  unobserved parents' configurations. Such \code{bn.fit} objects will be
  converted to $m$ code{grain} objects by replacing the missing conditional
  probability distributions with uniform distributions.

  Another solution to this problem is to fit another \code{bn.fit} with
  \code{method = "bayes"} and a low \code{iss} value, using the same data
  and network structure.

  Ordinal nodes will be treated as categorical by \code{as.grain},
  disregarding the ordering of the levels.

}
\value{

  An object of class \code{grain} (for \code{as.grain}), \code{bn.fit} (for
  \code{as.bn.fit}) or \code{bn} (for \code{as.bn}).

}
\examples{
\dontrun{
library(gRain)
a = bn.fit(hc(learning.test), learning.test)
b = as.grain(a)
c = as.bn.fit(b)}}
\author{Marco Scutari}
\keyword{utilities}
\keyword{interfaces to other packages}
