\name{bn.var}
\alias{bn.moments}
\alias{bn.var}
\title{ Structure variability of Bayesian networks }
\description{

  Measure the variability of the structure of a Bayesian network.

}
\usage{
# first and second moments' estimation
bn.moments(data, R = 200, m = nrow(data), algorithm,
  algorithm.args = list(), reduce = NULL, debug = FALSE)
# descriptive statistics
bn.var(x, method)
}
\arguments{
  \item{data}{a data frame containing the variables in the model.}
  \item{R}{a positive integer, the number of bootstrap replicates (in
    \code{bn.moments}) or the number of Monte Carlo samples (in
    \code{bn.var.test}).}
  \item{m}{a positive integer, the bootstrap sample size.}
  \item{algorithm}{a character string, the learning algorithm to be
    applied to the bootstrap replicates. Possible values are \code{gs},
    \code{iamb}, \code{fast.iamb}, \code{inter.iamb}, \code{mmpc},
    \code{hc}, \code{tabu}, \code{mmhc} and \code{rsmax2}. See
    \code{\link{bnlearn-package}} and the documentation of each algorithm
    for details.}
  \item{algorithm.args}{a list of extra arguments to be passed to
    the learning algorithm.}
  \item{x}{a covariance matrix or an object of class \code{mvber.moments}
    (the return value of the \code{bn.moments} function).}
  \item{method}{a character string, the label of the statistic. Possible
    values are \code{tvar} (\emph{total variance}), \code{gvar} 
    (\code{generalized variance}), \code{nvar} (\emph{Frobenius matrix
    norm}, which is equivalent to \emph{Nagao's test}) and \code{nvark}
    (another measure based on the \emph{Frobenius matrix norm}).}
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
}
}
\author{ Marco Scutari }
\keyword{multivariate}
