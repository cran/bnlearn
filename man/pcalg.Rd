\name{pcalg integration}
\alias{pcalg integration}
\alias{as.bn.pcAlgo}
\title{Import and export networks from the pcalg package}
\description{

  Convert \code{pcAlgo} objects to \code{bn} objects.

}
\usage{
\method{as.bn}{pcAlgo}(x, ..., check.cycles = TRUE)
}
\arguments{
  \item{x}{an object of class \code{pcAlgo}.}
  \item{\dots}{extra arguments from the generic method (currently ignored).}
  \item{check.cycles}{a boolean value. If \code{FALSE} the returned network will
    not be checked for cycles.}
}
\value{

  An object of class \code{bn}.

}
\author{Marco Scutari}
\keyword{utilities}
\keyword{interfaces to other packages}
