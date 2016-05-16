\name{score}
\alias{score}
\alias{logLik.bn}
\alias{AIC.bn}
\alias{BIC.bn}
\title{Score of the Bayesian network}
\description{

  Compute the score of the Bayesian network.

}
\usage{
score(x, data, type = NULL, ..., debug = FALSE)

\method{logLik}{bn}(object, data, ...)
\method{AIC}{bn}(object, data, ..., k = 1)
\method{BIC}{bn}(object, data, ...)
}
\arguments{
  \item{x, object}{an object of class \code{bn}.}
  \item{data}{a data frame containing the data the Bayesian network was learned
    from.}
  \item{type}{a character string, the label of a network score. If none is 
    specified, the default score is the \emph{Bayesian Information Criterion}
    for both discrete and continuous data sets. See 
    \code{\link{bnlearn-package}} for details.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
  \item{\dots}{extra arguments from the generic method (for the \code{AIC} and
    \code{logLik} functions, currently ignored) or additional tuning parameters
    (for the \code{score} function).}
  \item{k}{a numeric value, the penalty per parameter to be used; the default
    \code{k = 1} gives the expression used to compute the AIC in the context of
    scoring Bayesian networks.}
}
\details{

  Additional parameters of the \code{score} function:

  \itemize{

    \item \code{iss}: the imaginary sample size, used by the Bayesian Dirichlet
      equivalent score (both the \code{bde} and \code{mbde}) and the Bayesian
      Gaussian score (\code{bge}). It is also known as \dQuote{equivalent sample
      size}. The default value is equal to \code{10} for both the 
      \code{bde}/\code{mbde} scores and \code{bge}.
    \item \code{exp}: a list of indexes of experimental observations (those that
      have been artificially manipulated). Each element of the list must be
      named after one of the nodes, and must contain a numeric vector with 
      indexes of the observations whose value has been manipulated for that node.
    \item \code{k}: the penalty per parameter to be used by the AIC and BIC 
      scores. The default value is \code{1} for AIC and \code{log(nrow(data))/2}
      for BIC.
    \item \code{phi}: the prior phi matrix formula to use in the Bayesian
      Gaussian equivalent (\code{bge}) score. Possible values are 
      \code{heckerman} (default) and \code{bottcher} (the one used by default
      in the \pkg{deal} package.)
    \item \code{prior}: the prior distribution to be used with the various
      Bayesian Dirichlet scores (\code{bde}, \code{mbde}, \code{bds}) and the
      Bayesian Gaussian score (\code{bge}). Possible values are \code{uniform}
      (the default), \code{vsp} (the Bayesian variable selection prior, which
      puts a probability of inclusion on parents) and \code{cs} (the Castelo &
      Siebes prior, which puts an independent prior probability on each arc and
      direction).
    \item \code{beta}: the parameter associated with \code{prior}. If
      \code{prior} is \code{uniform}, \code{beta} is ignored. If \code{prior} is
      \code{vsp}, \code{beta} is the probability of inclusion of an additional
      parent (the default is \code{1/ncol(data)}). If \code{prior} is \code{cs},
      \code{beta} is a data frame with columns \code{from}, \code{to} and
      \code{prob} specifying the prior probability for a set of arcs. A uniform
      probability distribution is assumed for the remaining arcs.

  }

}
\note{

  AIC and BIC are computed as \code{logLik(x) - k * nparams(x)}, that is, the
  classic definition rescaled by -2. Therefore higher values are better, and
  for large sample sizes BIC converges to log(BDe).

  When using the Castelo & Siebes prior in structure learning, the prior
  probabilties associated with an arc are bound away from zero and one by
  shrinking them towards the uniform distribution as per Hausser and Strimmer
  (2009) with a lambda equal to \code{3 * sqrt(.Machine$double.eps)}. This 
  dramatically improves structure learning, which is less likely to get stuck
  when starting from an empty graph. As an alternative to prior probabilities,
  a blacklist can be used to prevent arcs from being included in the network,
  and a whitelist can be used to force the inclusion of particular arcs.
  \code{beta} is not modified when the prior is used from functions other than
  those implementing score-based and hybrid structure learning.

}
\references{

  Castelo R, Siebes A (2000). "Priors on Network Structures. Biasing the Search
      for Bayesian Networks". International Journal of Approximate Reasoning,
      \strong{24}(1), 39-57.

  Chickering DM (1995). "A Transformational Characterization of Equivalent
      Bayesian Network Structures". In "UAI '95: Proceedings of the Eleventh
      Annual Conference on Uncertainty in Artificial Intelligence", pp. 87-98.
      Morgan Kaufmann.

  Cooper GF, Yoo C (1999). "Causal Discovery from a Mixture of Experimental
      and Observational Data". In "UAI '99: Proceedings of the Fifteenth
      Annual Conference on Uncertainty in Artificial Intelligence", pp. 116-125.
      Morgann Kaufmann.

  Geiger D, Heckerman D (1994). "Learning Gaussian Networks". In "UAI '94: 
      Proceedings of the Tenth Annual Conference on Uncertainty in Artificial
      Intelligence", pp. 235-243. Morgann Kaufmann.
      Available as Technical Report MSR-TR-94-10.

  Hausser J, Strimmer K (2009). "Entropy inference and the James-Stein
      estimator, with application to nonlinear gene association networks".
      \emph{Statistical Applications in Genetics and Molecular Biology},
      \strong{10}, 1469-1484.

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

## let's see score equivalence in action!
res2 = set.arc(gs(learning.test), "B", "A")
score(res2, learning.test, type = "bde")

## K2 score on the other hand is not score equivalent.
score(res, learning.test, type = "k2")
score(res2, learning.test, type = "k2")

## BDe with a prior.
beta = data.frame(from = c("A", "D"), to = c("B", "F"),
         prob = c(0.2, 0.5), stringsAsFactors = FALSE)
score(res, learning.test, type = "bde", prior = "cs", beta = beta)

## equivalent to logLik(res, learning.test)
score(res, learning.test, type = "loglik")

## equivalent to AIC(res, learning.test)
score(res, learning.test, type = "aic")
}
\seealso{\code{\link{choose.direction}}, \code{\link{arc.strength}},
  \code{\link{alpha.star}}.}
\author{ Marco Scutari }
\keyword{htest}
