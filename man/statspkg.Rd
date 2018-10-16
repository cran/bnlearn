\name{lm integration}
\alias{lm integration}
\alias{as.lm}
\alias{as.lm.bn}
\alias{as.lm.bn.fit}
\alias{as.lm.bn.fit.gnode}
\title{Produce lm objects from Bayesian networks}
\description{

  Take a \code{bn} object or \code{bn.fit} object encoding a Gaussian network
  and refit all the local distributions using \code{lm()}. This makes it
  possible to use all the functions provided by R for \code{lm} objects 
  (\code{summary}, \code{anova}, etc.) to investigate the network.

}
\usage{
\method{as.lm}{bn}(x, data, ...)
\method{as.lm}{bn.fit}(x, data, ...)
\method{as.lm}{bn.fit.gnode}(x, data, ...)
}
\arguments{
  \item{x}{an object of class \code{bn}, \code{bn.fit} or \code{bn.fit.gnode}.}
  \item{data}{a data frame containing the variables in the model.}
  \item{\dots}{additional arguments, currently ignored.}
}
\value{

  If \code{x} is an object of class \code{bn} or \code{bn.fit}, \code{as.lm()}
  returns a list of \code{lm} objects, one for each node in \code{x}. If
  \code{x} is an object of class \code{bn} or \code{bn.fit.gnode},
  \code{as.lm()} returns a single \code{lm} object.

}
\examples{
dag = hc(gaussian.test)
fitted = bn.fit(dag, gaussian.test)
as.lm(dag, gaussian.test)
as.lm(fitted, gaussian.test)
as.lm(fitted$F, gaussian.test)
}
\author{Marco}
\keyword{interfaces to other packages}
