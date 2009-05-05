\name{score}
\alias{score}
\alias{logLik.bn}
\alias{AIC.bn}
\title{ Score of the Bayesian network }
\description{

  Compute the score of the Bayesian network.

}
\usage{
  score(x, data, type = NULL, ..., debug = FALSE)

  \method{logLik}{bn}(object, data, ...)
  \method{AIC}{bn}(object, data, ..., k = 1)
}
\arguments{
  \item{x}{an object of class \code{bn}.}
  \item{object}{an object of class \code{bn}.}
  \item{data}{a data frame, containing the data the Bayesian network
       was learned from.}
  \item{type}{a character string, the label of a network score. 
      If none is specified, the default score is the 
      \emph{Akaike Information Criterion} for discrete data sets 
      and the \emph{Bayesian Gaussian posterior density} for 
      continuous ones. See \code{\link{bnlearn-package}}
      for details.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output
       is printed; otherwise the function is completely silent.}
  \item{\dots}{extra arguments from the generic method (for the \code{AIC}
       and \code{logLik} functions, currently ignored) or additional
       tuning parameters (for the \code{score} function).}
  \item{k}{a numeric value, the penalty per parameter to be used; the
       default \code{k = 1} gives the expression used to compute the
       AIC in the context of scoring Bayesian networks.}
}
\details{

  Additional parameters of the \code{score} function:

  \itemize{

    \item \code{iss}: the imaginary sample size, used by the Bayesian
        Dirichlet equivalent score and the Bayesian Gaussian posterior
        density. The default value is equal to the number of cells of the
        joint contingency table (for compatibility with the \pkg{deal}
        package) for the \code{bde} score, or to the number of
        independent parameters for the \code{bge} score.

    \item \code{k}: the penalty per parameter to be used by the AIC and
        BIC scores. The default value is \code{1} for AIC and
        \code{log(nrow(data))/2} for BIC.

    \item \code{phi}: the prior phi matrix formula to use in the
        Bayesian Gaussian equivalent (\code{bge}) score. Possible
        values are \code{heckerman} (default) and \code{bottcher}
        (the one used by default in the \pkg{deal} package.)

  }

}
\note{

  Execution time scales linearly with the sample size.

}
\references{

  Chickering DM (1995). "A Transformational Characterization of Equivalent 
      Bayesian Network Structures". In "UAI '95: Proceedings of the Eleventh 
      Annual Conference on Uncertainty in Artificial Intelligence", pp. 87-98.
      Morgan Kaufmann.

  Geiger D, Heckerman D (1994). "Learning Gaussian Networks". \emph{Technical 
      report}, Microsoft Research. Available as Technical Report MSR-TR-94-10.

  Heckerman D, Geiger D, Chickering DM (1995). "Learning Bayesian Networks: The
      Combination of Knowledge and Statistical Data". \emph{Machine Learning}, 
      \strong{20}(3), 197-243. Available as Technical Report MSR-TR-94-09.

}
\value{

  A numeric value, the score of the Bayesian network.

}
\examples{
data(learning.test)
res = set.arc(gs(learning.test), "A", "B")
score(res, learning.test, type = "bde")
# [1] -24002.36
## let's see score equivalence in action!
res2 = set.arc(gs(learning.test), "B", "A")
score(res2, learning.test, type = "bde")
# [1] -24002.36

## k2 score on the other hand is not score equivalent.
score(res, learning.test, type = "k2")
# [1] -23958.70
score(res2, learning.test, type = "k2")
# [1] -23957.68

## equivalent to logLik(res, learning.test)
score(res, learning.test, type = "loglik")
# [1] -23832.13

## equivalent to AIC(res, learning.test)
score(res, learning.test, type = "aic")
# [1] -23873.13
}
\seealso{\code{\link{choose.direction}}, \code{\link{arc.strength}}. }
\author{ Marco Scutari }
\keyword{htest}
