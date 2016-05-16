\name{bn.fit utilities}
\alias{bn.fit utilities}
\alias{sigma}
\alias{fitted.bn.fit}
\alias{coef.bn.fit}
\alias{residuals.bn.fit}
\alias{sigma.bn.fit}
\alias{predict.bn.fit}
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
\method{predict}{bn.fit}(object, node, data, method = "parents", ..., prob = FALSE,
  debug = FALSE)
\method{logLik}{bn.fit}(object, data, nodes, by.sample = FALSE, ...)
\method{AIC}{bn.fit}(object, data, ..., k = 1)
\method{BIC}{bn.fit}(object, data, ...)

## methods available for "bn.fit.dnode"
\method{coef}{bn.fit.dnode}(object, ...)

## methods available for "bn.fit.onode"
\method{coef}{bn.fit.onode}(object, ...)

## methods available for "bn.fit.gnode"
\method{fitted}{bn.fit.gnode}(object, ...)
\method{coef}{bn.fit.gnode}(object, ...)
\method{residuals}{bn.fit.gnode}(object, ...)
\method{sigma}{bn.fit.gnode}(object, ...)

## methods available for "bn.fit.cgnode"
\method{fitted}{bn.fit.cgnode}(object, ...)
\method{coef}{bn.fit.cgnode}(object, ...)
\method{residuals}{bn.fit.cgnode}(object, ...)
\method{sigma}{bn.fit.cgnode}(object, ...)
}
\arguments{
  \item{object}{an object of class \code{bn.fit}, \code{bn.fit.dnode},
    \code{bn.fit.gnode}, \code{bn.fit.cgnode} or \code{bn.fit.onode}.}
  \item{node}{a character string, the label of a node.}
  \item{nodes}{a vector of character strings, the label of a nodes whose
    log-likelihood components are to be computed.}
  \item{data}{a data frame containing the variables in the model.}
  \item{method}{a character string, the method used to estimate predictions.
    See below.}
  \item{\dots}{additional arguments. See below.}
  \item{k}{a numeric value, the penalty per parameter to be used; the default
    \code{k = 1} gives the expression used to compute AIC.}
  \item{by.sample}{a boolean value. If \code{TRUE}, \code{logLik} returns a
    vector containing the the log-likelihood of each observations in the
    sample. If \code{FALSE}, \code{logLik} returns a single value, the 
    likelihood of the whole sample.}
  \item{prob}{a boolean value. If \code{TRUE} and \code{object} is a discrete
    network, the probabilities used for prediction are attached to the
    predicted values as an attribute called \code{prob}.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}

}
\details{

  \code{coef} (and its alias \code{coefficients}) extracts model coefficients
  (which are conditional probabilities for discrete nodes and linear regression
  coefficients for Gaussian and conditional Gaussian nodes).

  \code{residuals} (and its alias \code{resid}) extracts model residuals and
  \code{fitted} (and its alias \code{fitted.values}) extracts fitted values
  from Gaussian and conditional Gaussian nodes. If the \code{bn.fit} object
  does not include the residuals or the fitted values for the node of interest
  both functions return \code{NULL}.

  \code{sigma} extracts the standard deviations of the residuals from Gaussian
  and conditional Gaussian networks and nodes.

  \code{logLik} returns the log-likelihood for the observations in \code{data}.

  \code{predict} returns the predicted values for \code{node} given the data
  specified by \code{data} and the fitted network. Depending on the value of
  \code{method}, the predicted values are computed as follows.
  \itemize{
    \item \code{parents}: the predicted values are computed by plugging in
      the new values for the parents of \code{node} in the local probability
      distribution of \code{node} extracted from \code{fitted}.
    \item \code{bayes-lw}: the predicted values are computed by averaging
      likelihood weighting simulations performed using all the available nodes
      as evidence (obviously, with the exception of the node whose values we
      are predicting). The number of random samples which are averaged for each
      new observation is controlled by the \code{n} optional argument; the
      default is \code{500}. If the variable being predicted is discrete, the
      predicted level is that with the highest  conditional probability. If the
      variable is continuous, the predicted value is the expected value of the
      conditional distribution.
  }

}
\note{

  Ties in prediction are broken using \emph{Bayesian tie breaking}, i.e.
  sampling at random from the tied values. Therefore, setting the random seed
  is required to get reproducible results.

  \code{predict} accepts either a \code{bn} or a \code{bn.fit} object as its
  first argument. For the former, the parameters of the network are fitted on
  \code{data}, that is, the observations whose class labels the function is
  trying to predict.

}
\value{

  \code{predict} returns a numeric vector (for Gaussian and conditional
  Gaussian nodes), a factor (for categorical nodes) or an ordered factor (for
  ordinal nodes). If \code{prob = TRUE} and the network is discrete, the
  probabilities used for prediction are attached to the predicted values as
  an attribute called \code{prob}.

  \code{logLik} returns a numeric vector or a single numeric value, depending
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
}
\author{ Marco Scutari }
\keyword{utilities}
\keyword{graphs}
