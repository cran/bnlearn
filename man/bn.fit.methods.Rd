\name{bn.fit utilities}
\alias{bn.fit utilities}
\alias{fitted.bn.fit}
\alias{coef.bn.fit}
\alias{residuals.bn.fit}
\alias{predict.bn.fit}
\alias{logLik.bn.fit}
\alias{AIC.bn.fit}
\alias{BIC.bn.fit}
\alias{fitted.bn.fit.gnode}
\alias{coef.bn.fit.gnode}
\alias{residuals.bn.fit.gnode}
\alias{predict.bn.fit.gnode}
\alias{fitted.bn.fit.dnode}
\alias{coef.bn.fit.dnode}
\alias{residuals.bn.fit.dnode}
\alias{predict.bn.fit.dnode}
\alias{coef.bn.fit.onode}
\alias{predict.bn.fit.onode}

\title{ Utilities to manipulate fitted Bayesian networks }
\description{

 Assign, extract or compute various quantities of interest from an object
 of class \code{bn.fit}, \code{bn.fit.dnode}, \code{bn.fit.gnode} or
 \code{bn.fit.onode}.

}
\usage{
## methods available for "bn.fit"
\method{fitted}{bn.fit}(object, ...)
\method{coef}{bn.fit}(object, ...)
\method{residuals}{bn.fit}(object, ...)
\method{predict}{bn.fit}(object, node, data, ..., debug = FALSE)
\method{logLik}{bn.fit}(object, data, nodes, by.sample = FALSE, ...)
\method{AIC}{bn.fit}(object, data, ..., k = 1)
\method{BIC}{bn.fit}(object, data, ...)

## methods available for "bn.fit.dnode"
\method{coef}{bn.fit.dnode}(object, ...)
\method{predict}{bn.fit.dnode}(object, data, ..., debug = FALSE)

## methods available for "bn.fit.onode"
\method{coef}{bn.fit.onode}(object, ...)
\method{predict}{bn.fit.onode}(object, data, ..., debug = FALSE)

## methods available for "bn.fit.gnode"
\method{fitted}{bn.fit.gnode}(object, ...)
\method{coef}{bn.fit.gnode}(object, ...)
\method{residuals}{bn.fit.gnode}(object, ...)
\method{predict}{bn.fit.gnode}(object, data, ..., debug = FALSE)
}
\arguments{
  \item{object}{an object of class \code{bn.fit}, \code{bn.fit.dnode}
    or \code{bn.fit.gnode}.}
  \item{node}{a character string, the label of a node.}
  \item{nodes}{a vector of character strings, the label of a nodes whose
        loglikelihood components are to be computed.}
  \item{data}{a data frame containing the variables in the model.}
  \item{\dots}{additional arguments (currently ignored).}
  \item{k}{a numeric value, the penalty per parameter to be used; the
       default \code{k = 1} gives the expression used to compute AIC.}
  \item{by.sample}{a boolean value. If \code{TRUE}, \code{logLik} returns
       a vector containing the the log-likelihood of each observations in
       the sample. If \code{FALSE}, \code{logLik} returns a single value,
       the likelihood of the whole sample.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output
       is printed; otherwise the function is completely silent.}

}
\details{

  \code{coef} (and its alias \code{coefficients}) extracts model
  coefficients (which are conditional probabilities in discrete
  networks and linear regression coefficients in Gaussian networks).

  \code{residuals} (and its alias \code{resid}) extracts model
  residuals and \code{fitted} (and its alias \code{fitted.values})
  extracts fitted values from fitted Gaussian networks.
  If the \code{bn.fit} object does not include the residuals or the
  fitted values (for the nodes of interest, in the case of
  \code{bn.fit.gnode} objects), both functions return \code{NULL}.

  \code{predict} returns the predicted values for \code{node} for
  the data specified by \code{data}.

}
\note{

  Ties in prediction are broken using \emph{Bayesian tie breaking}, i.e.
  sampling at random from the tied values. Therefore, setting the random
  seed is required to get reproducible results.

}
\value{

  \code{predict} returns a numeric vector (for Gaussian networks) or a factor
  (for discrete networks).

  \code{logLik} returns a numeric vector or a single numeric value, depending
  on the value of \code{by.sample}. \code{AIC} and \code{BIC} always return
  a single numeric value.

  All the other functions return a list with an element for each node in
  the network (if \code{object} has class \code{bn.fit}) or a numeric
  vector (if \code{object} has class \code{bn.fit.dnode} or \code{bn.fit.gnode}).

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
}
\author{ Marco Scutari }
\keyword{utilities}
\keyword{graphs}
