\name{bn.cv}
\alias{bn.cv}
\title{ Cross-validation for Bayesian networks }
\description{

  Perform a k-fold cross-validation for a learning algorithm
  or a fixed network structure.

}
\usage{
bn.cv(data, bn, loss = NULL, k = 10, algorithm.args = list(),
  loss.args = list(), fit = "mle", fit.args = list(),
  cluster = NULL, debug = FALSE)
}
\arguments{
  \item{data}{a data frame containing the variables in the model.}
  \item{bn}{either a character string (the label of the learning
    algorithm to be applied to the training data in each iteration)
    or an object of class \code{bn} (a fixed network structure).}
  \item{loss}{a character string, the label of a loss function.
    If none is specified, the default loss function is the
    \emph{Log-Likelihood Loss} for both discrete and continuous
    data sets. See below for additional details.}
  \item{k}{a positive integer number, the number of groups into
    which the data will be split.}
  \item{algorithm.args}{a list of extra arguments to be passed to
    the learning algorithm.}
  \item{loss.args}{a list of extra arguments to be passed to
    the loss function specified by \code{loss}.}
  \item{fit}{a character string, the label of the method used to fit the
    parameters of the newtork. See \code{\link{bn.fit}} for details.}
  \item{fit.args}{additional arguments for the parameter estimation prcoedure,
    see again \code{\link{bn.fit}} for details..}
  \item{cluster}{an optional cluster object from package \pkg{snow}.
      See \code{\link{snow integration}} for details and a simple
      example.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output
    is printed; otherwise the function is completely silent.}
}
\details{

  The following loss functions are implemented:

  \itemize{

    \item \emph{Log-Likelihood Loss} (\code{logl}): also known as \emph{negative
      entropy} or \emph{negentropy}, it's the negated expected log-likelihood
      of the test set for the Bayesian network fitted from the training set.

    \item \emph{Gaussian Log-Likelihood Loss} (\code{logl-g}): the negated expected
      log-likelihood for Gaussian Bayesian networks.

    \item \emph{Classification Error} (\code{pred}): the \emph{prediction error}
      for a single node (specified by the \code{target} parameter in \code{loss.args})
      in a discrete network.

  }

}
\value{

  An object of class \code{bn.kcv}.

}
\author{ Marco Scutari }
\seealso{\code{\link{bn.boot}}, \code{\link{rbn}}, \code{\link{bn.kcv-class}}.}
\references{

  Koller D, Friedman N (2009). \emph{Probabilistic Graphical Models:
      Principles and Techniques}. MIT Press.

}
\examples{
bn.cv(learning.test, 'hc', loss = "pred", loss.args = list(target = "F"))
#
#  k-fold cross-validation for Bayesian networks
#
#  target learning algorithm:             Hill-Climbing
#  number of subsets:                     10
#  loss function:                         Classification Error
#  expected loss:                         0.509
#
bn.cv(gaussian.test, 'mmhc')
#
#  k-fold cross-validation for Bayesian networks
#
#  target learning algorithm:             Max-Min Hill Climbing
#  number of subsets:                     10
#  loss function:                         Log-Likelihood Loss (Gaussian)
#  expected loss:                         10.63062
#
}
\keyword{nonparametric}
\keyword{multivariate}
