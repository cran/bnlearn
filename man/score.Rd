\name{score}
\alias{score}
\alias{score,bn-method}
\alias{score,bn.naive-method}
\alias{score,bn.tan-method}
\alias{logLik.bn}
\alias{AIC.bn}
\alias{BIC.bn}
\title{Score of the Bayesian network}
\description{

  Compute the score of the Bayesian network.

}
\usage{
\S4method{score}{bn}(x, data, type = NULL, ..., by.node = FALSE, debug = FALSE)
\S4method{score}{bn.naive}(x, data, type = NULL, ..., by.node = FALSE, debug = FALSE)
\S4method{score}{bn.tan}(x, data, type = NULL, ..., by.node = FALSE, debug = FALSE)

\method{logLik}{bn}(object, data, ...)
\method{AIC}{bn}(object, data, ..., k = 1)
\method{BIC}{bn}(object, data, ...)
}
\arguments{
  \item{x, object}{an object of class \code{bn}.}
  \item{data}{a data frame containing the data the Bayesian network that will
    be used to compute the score.}
  \item{type}{a character string, the label of a network score. If none is
    specified, the default score is the \emph{Bayesian Information Criterion}
    for both discrete and continuous data sets. See \code{\link{network scores}}
    for details.}
  \item{by.node}{a boolean value. If \code{TRUE} and the score is decomposable,
    the function returns the score terms corresponding to each node; otherwise
    it returns their sum (the overall score of \code{x}).}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
  \item{\dots}{extra arguments from the generic method (for the \code{AIC} and
    \code{logLik} functions, currently ignored) or additional tuning parameters
    (for the \code{score} function).}
  \item{k}{a numeric value, the penalty coefficient to be used; the default
    \code{k = 1} gives the expression used to compute the AIC in the context of
    scoring Bayesian networks.}
}
\details{

  Additional arguments of the \code{score()} function:

  \itemize{

    \item \code{iss}: the imaginary sample size used by the Bayesian Dirichlet
      scores (\code{bde}, \code{mbde}, \code{bds}, \code{bdj}). It is also known
      as \dQuote{equivalent sample size}. The default value is equal to \code{1}.
    \item \code{iss.mu}: the imaginary sample size for the normal component of
      the normal-Wishart prior in the Bayesian Gaussian score (\code{bge}). The
      default value is \code{1}.
    \item \code{iss.w}: the imaginary sample size for the Wishart component of
      the normal-Wishart prior in the Bayesian Gaussian score (\code{bge}). The
      default value is \code{ncol(data) + 2}.
    \item \code{nu}: the mean vector of the normal component of the
      normal-Wishart prior in the Bayesian Gaussian score (\code{bge}). The
      default value is equal to \code{colMeans(data)}.
    \item \code{l}: the number of scores to average in the locally averaged
      Bayesian Dirichlet score (\code{bdla}). The default value is \code{5}.
    \item \code{exp}: a list of indexes of experimental observations (those that
      have been artificially manipulated). Each element of the list must be
      named after one of the nodes, and must contain a numeric vector with
      indexes of the observations whose value has been manipulated for that node.
    \item \code{k}: the penalty coefficient to be used by the AIC and BIC
      scores. The default value is \code{1} for AIC and \code{log(nrow(data))/2}
      for BIC.
    \item \code{prior}: the prior distribution to be used with the various
	  Bayesian Dirichlet scores (\code{bde}, \code{mbde}, \code{bds}, \code{bdj},
	  \code{bdla}) and the Bayesian Gaussian score (\code{bge}). Possible
	  values are \code{uniform} (the default), \code{vsp} (the Bayesian
	  variable selection prior, which puts a probability of inclusion on
	  parents), \code{marginal} (an independent marginal uniform for each arc)
	  and \code{cs} (the Castelo & Siebes prior, which puts an independent
	  prior probability on each arc and direction).
    \item \code{beta}: the parameter associated with \code{prior}.

      \itemize{

        \item If \code{prior} is \code{uniform}, \code{beta} is ignored.
        \item If \code{prior} is \code{vsp}, \code{beta} is the probability of
          inclusion of an additional parent. The default is \code{1/ncol(data)}.
        \item If \code{prior} is \code{marginal}, \code{beta} is the probability
          of inclusion of an arc. Each direction has a probability of inclusion
          of \code{beta / 2} and the probability that the arc is not included is
          therefore  \code{1 - beta}. The default value is \code{0.5}, so that
          arc inclusion and arc exclusion have the same probability.
        \item If \code{prior} is \code{cs}, \code{beta} is a data frame with
          columns \code{from}, \code{to} and \code{prob} specifying the prior
          probability for a set of arcs. A uniform probability distribution is
          assumed for the remaining arcs.

      }
    \item \code{newdata}: the test set whose predictive likelihood will be
      computed by \code{pred-loglik}, \code{pred-loglik-g} or
      \code{pred-loglik-cg}. It should be a data frame with the same variables
      as \code{data}.

  }

}
\note{

  AIC and BIC are computed as \code{logLik(x) - k * nparams(x)}, that is, the
  classic definition rescaled by -2. Therefore higher values are better, and
  for large sample sizes BIC converges to log(BDe).

  When using the Castelo & Siebes prior in structure learning, the prior
  probabilities associated with an arc are bound away from zero and one by
  shrinking them towards the uniform distribution as per Hausser and Strimmer
  (2009) with a lambda equal to \code{3 * sqrt(.Machine$double.eps)}. This
  dramatically improves structure learning, which is less likely to get stuck
  when starting from an empty graph. As an alternative to prior probabilities,
  a blacklist can be used to prevent arcs from being included in the network,
  and a whitelist can be used to force the inclusion of particular arcs.
  \code{beta} is not modified when the prior is used from functions other than
  those implementing score-based and hybrid structure learning.

}
\value{

  For \code{score()} with \code{by.node = TRUE}, a vector of numeric values, the
  individual node contributions to the score of the Bayesian network.
  Otherwise, a single numeric value, the score of the Bayesian network.

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
\author{Marco Scutari}
\keyword{network scores}
