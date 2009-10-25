\name{bn.var}
\alias{bn.moments}
\alias{bn.var}
\alias{bn.var.test}
\title{ Structure Variability of Bayesian networks }
\description{

  Measure the variability of the structure of a Bayesian network.

}
\usage{
# first and second moments' estimation
bn.moments(data, R = 200, m = nrow(data), algorithm,
  algorithm.args = list(), reduce = NULL, debug = FALSE)
# descriptive statistics
bn.var(x, method)
# Monte Carlo test for entropy
bn.var.test(x, method, R, B, debug = FALSE)
}
\arguments{
  \item{data}{a data frame, containing the variables in the model.}
  \item{R}{a positive integer, the number of bootstrap replicates (in
    \code{bn.moments}) or the number of Monte Carlo samples (in
    \code{bn.var.test}).}
  \item{m, B}{a positive integer, the size of each bootstrap (in
    \code{bn.moments}) or Monte Carlo (in \code{bn.var.test}) replicate.}
  \item{algorithm}{a character string, the learning algorithm to be
    applied to the bootstrap replicates. Possible values are \code{gs},
    \code{iamb}, \code{fast.iamb}, \code{inter.iamb}, \code{mmpc}
    and \code{hc}. See \code{\link{bnlearn-package}} and the
    documentation of each algorithm for details.}
  \item{algorithm.args}{a list of extra arguments to be passed to
    the learning algorithm.}
  \item{x}{a covariance matrix or an object of class \code{mvber.moments}
    (the return value of the \code{bn.moments} function).}
  \item{method}{a character string, the label of the statistic used
    in \code{bn.var} or \code{bn.var.test}. Possible values are
    \code{tvar} (\emph{total variance}), \code{gvar} (\code{generalized
    variance}), \code{nvar} (\emph{Frobenius matrix norm}, which is
    equivalent to \emph{Nagao's test}).}
  \item{reduce}{a character string, either \code{first} or \code{second}.
    If \code{first} all the arcs with first moment equal to zero are
    dropped; if if \code{second} all the arcs with zero variance
    are dropped.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output
       is printed; otherwise the function is completely silent.}
}
\note{

  These functions are experimental implementations of techniques still
  in development; their form (name, parameters, etc.) will likely
  change without notice in the future.

}
\value{

  \code{bn.moments} returns an object of class \code{mvber.moments}.

  \code{bn.var} returns a vector of two elements, the observed value of
  the statistic (named \code{statistic}) and its normalized equivalent
  (named \code{normalized}).

  \code{bn.var.test} returns an object of class \code{htest}.

}
\references{

  Scutari M (2009). "Structure Variability in Bayesian Networks".
      \emph{ArXiv Statistics - Methodology e-prints}.
      http://arxiv.org/abs/0909.1685.

}
\examples{
\dontrun{
z =  bn.moments(learning.test, algorithm = "gs", R = 100)
bn.var(z, method = "tvar")
#  statistic normalized
#    1.29060    0.34416
bn.var.test(z, method = "nvar")
#
# 	squared Frobenius norm
#
# data:  covariance matrix
# nvar = 0.5471, B = 5000, R = 100, p-value < 2.2e-16
# alternative hypothesis: true value is greater than 0
}
}
\author{ Marco Scutari }
\keyword{utilities}
