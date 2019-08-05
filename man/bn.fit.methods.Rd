\name{bn.fit utilities}
\alias{bn.fit utilities}
\alias{sigma}
\alias{fitted.bn.fit}
\alias{coef.bn.fit}
\alias{residuals.bn.fit}
\alias{sigma.bn.fit}
\alias{logLik.bn.fit}
\alias{AIC.bn.fit}
\alias{BIC.bn.fit}
\alias{fitted.bn.fit.gnode}
\alias{coef.bn.fit.gnode}
\alias{residuals.bn.fit.gnode}
\alias{sigma.bn.fit.gnode}
\alias{fitted.bn.fit.dnode}
\alias{coef.bn.fit.dnode}
\alias{residuals.bn.fit.dnode}
\alias{coef.bn.fit.onode}
\alias{fitted.bn.fit.cgnode}
\alias{coef.bn.fit.cgnode}
\alias{residuals.bn.fit.cgnode}
\alias{sigma.bn.fit.cgnode}
\title{Utilities to manipulate fitted Bayesian networks}
\description{

 Assign, extract or compute various quantities of interest from an object of
 class \code{bn.fit}, \code{bn.fit.dnode}, \code{bn.fit.gnode},
 \code{bn.fit.cgnode} or \code{bn.fit.onode}.

}
\usage{
## methods available for "bn.fit"
\method{fitted}{bn.fit}(object, ...)
\method{coef}{bn.fit}(object, ...)
\method{residuals}{bn.fit}(object, ...)
\method{sigma}{bn.fit}(object, ...)
\method{logLik}{bn.fit}(object, data, nodes, by.sample = FALSE, ...)
\method{AIC}{bn.fit}(object, data, ..., k = 1)
\method{BIC}{bn.fit}(object, data, ...)

## methods available for "bn.fit.dnode"
\method{coef}{bn.fit.dnode}(object, for.parents, ...)

## methods available for "bn.fit.onode"
\method{coef}{bn.fit.onode}(object, for.parents, ...)

## methods available for "bn.fit.gnode"
\method{fitted}{bn.fit.gnode}(object, ...)
\method{coef}{bn.fit.gnode}(object, ...)
\method{residuals}{bn.fit.gnode}(object, ...)
\method{sigma}{bn.fit.gnode}(object, ...)

## methods available for "bn.fit.cgnode"
\method{fitted}{bn.fit.cgnode}(object, ...)
\method{coef}{bn.fit.cgnode}(object, for.parents, ...)
\method{residuals}{bn.fit.cgnode}(object,  ...)
\method{sigma}{bn.fit.cgnode}(object, for.parents, ...)
}
\arguments{
  \item{object}{an object of class \code{bn.fit}, \code{bn.fit.dnode},
    \code{bn.fit.gnode}, \code{bn.fit.cgnode} or \code{bn.fit.onode}.}
  \item{nodes}{a vector of character strings, the label of a nodes whose
    log-likelihood components are to be computed.}
  \item{data}{a data frame containing the variables in the model.}
  \item{\dots}{additional arguments, currently ignored.}
  \item{k}{a numeric value, the penalty coefficient to be used; the default
    \code{k = 1} gives the expression used to compute AIC.}
  \item{by.sample}{a boolean value. If \code{TRUE}, \code{logLik()} returns
    a vector containing the the log-likelihood of each observations in the
    sample. If \code{FALSE}, \code{logLik()} returns a single value, the
    likelihood of the whole sample.}
  \item{for.parents}{a named list in which each element contains }

}
\details{

  \code{coef()} (and its alias \code{coefficients()}) extracts model coefficients
  (which are conditional probabilities for discrete nodes and linear regression
  coefficients for Gaussian and conditional Gaussian nodes).

  \code{residuals()} (and its alias \code{resid()}) extracts model residuals and
  \code{fitted()} (and its alias \code{fitted.values()}) extracts fitted values
  from Gaussian and conditional Gaussian nodes. If the \code{bn.fit} object
  does not include the residuals or the fitted values for the node of interest
  both functions return \code{NULL}.

  \code{sigma()} extracts the standard deviations of the residuals from Gaussian
  and conditional Gaussian networks and nodes.

  \code{logLik()} returns the log-likelihood for the observations in \code{data}.

  The \code{for.parents} argument in the methods for \code{coef()} and
  \code{sigma()} can be used to have both functions return the parameters
  associated with a specific configuration of the discrete parents of a node.
  If \code{for.parents} is not specified, all relevant parameters are returned.

}
\value{

  \code{logLik()} returns a numeric vector or a single numeric value, depending
  on the value of \code{by.sample}. \code{AIC} and \code{BIC} always return a
  single numeric value.

  All the other functions return a list with an element for each node in the
  network (if \code{object} has class \code{bn.fit}) or a numeric vector or
  matrix (if \code{object} has class \code{bn.fit.dnode}, \code{bn.fit.gnode},
  \code{bn.fit.cgnode} or \code{bn.fit.onode}).

}
\seealso{
  \code{\link{bn.fit}}, \code{\link{bn.fit-class}}.
}
\examples{
data(gaussian.test)
res = hc(gaussian.test)
fitted = bn.fit(res, gaussian.test)
coefficients(fitted)
coefficients(fitted$C)
str(residuals(fitted))

data(learning.test)
res2 = hc(learning.test)
fitted2 = bn.fit(res2, learning.test)
coefficients(fitted2$E)
coefficients(fitted2$E, for.parents = list(F = "a", B = "b"))
}
\author{Marco Scutari}
\keyword{network scores}
\keyword{convenience functions}
