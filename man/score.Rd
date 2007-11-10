\name{score}
\alias{score}
\alias{logLik.bn}
\alias{AIC.bn}
\title{ Score of the Bayesian network }
\description{

  Compute the score of the Bayesian network.

}
\usage{

  score(x, data, type = "dir", debug = FALSE)

  \method{logLik}{bn}(object, data, ...)
  \method{AIC}{bn}(object, data, ..., k = 1)

}
\arguments{
  \item{x}{an object of class "bn".}
  \item{object}{an object of class "bn".}
  \item{data}{a data frame, containing the data the Bayesian network
       was learned from.}
  \item{type}{a character string, the label of the score. Valid values 
       are \code{lik} (likelihood), \code{loglik} (log-likelihood), 
       \code{aic} (Akaike Information Criterion aka AIC), \code{bic} 
       (Bayesian Information Criterion aka BIC) and \code{dir} 
       (Bayesian Dirichlet aka BDe).}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output 
       is printed; otherwise the function is completely silent.}
  \item{...}{extra arguments from the generic method (currently ignored).}
  \item{k}{a numeric value, the penalty per parameter to be used; the 
          default \code{k = 1} gives the expression used to compute the
          AIC in the context of scoring Bayesian networks.}
}
\note{

  Only discrete Bayesian networks are supported.

  The execution time scales linearly with the sample size.

}
\references{

  D. M. Chickering. A Transformational Characterization of Equivalent Bayesian
    Network Structures. In Proceedins of 11th Conference on Uncertainty in
    Artificial Intelligence, pages 87-98. Morgan Kaufmann Publishers Inc., 1995.

}
\value{

  A numeric value, the score of the Bayesian network.

}
\examples{

data(learning.test)
res = set.arc(gs(learning.test), "A", "B")
score(res, learning.test, type = "dir")
# [1] -25005.24
# let's see score equivalence in action!
res = set.arc(gs(learning.test), "B", "A")
score(res, learning.test, type = "dir")
# [1] -25005.24

# equivalent to logLik(res, learning.test)
score(res, learning.test, type = "loglik")
# [1] -129.9684

# equivalent to AIC(res, learning.test)
score(res, learning.test, type = "aic")
# [1] -170.9684

}
\author{ Marco Scutari }
\keyword{htest}
