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
score(x, data, type = NULL, ..., by.node = FALSE, debug = FALSE)

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
    for both discrete and continuous data sets. See
    \code{\link{bnlearn-package}} for details.}
  \item{by.node}{a boolean value. If \code{TRUE} and the score is decomposable,
    the function returns the individual node contributions to the score;
    otherwise it returns their sum (the overall score of \code{x}).}
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

  Additional arguments of the \code{score} function:

  \itemize{

    \item \code{iss}: the imaginary sample size, used by the Bayesian Dirichlet
      scores (\code{bde}, \code{mbde}, \code{bds}, \code{bdj}) and the Bayesian
      Gaussian score (\code{bge}). It is also known as \dQuote{equivalent sample
      size}. The default value is equal to \code{1} for the Dirichlet scores and
      \code{10} for \code{bge}.
    \item \code{exp}: a list of indexes of experimental observations (those that
      have been artificially manipulated). Each element of the list must be
      named after one of the nodes, and must contain a numeric vector with
      indexes of the observations whose value has been manipulated for that node.
    \item \code{k}: the penalty coefficient to be used by the AIC and BIC
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
      puts a probability of inclusion on parents), \code{marginal} (an
      independent marginal uniform for each arc) and \code{cs} (the Castelo &
      Siebes prior, which puts an independent prior probability on each arc
      and direction).
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
      for Bayesian Networks". \emph{International Journal of Approximate
      Reasoning}, \strong{24}(1), 39-57.

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

  Suzuki J (2016). "A Theoretical Analysis of the BDeu Scores in Bayesian
      Network Structure Learning". \emph{Behaviormetrika}, \strong{44}(1),
      97-116. 

  Scutari M (2016). "An Empirical-Bayes Score for Discrete Bayesian Networks".
      \emph{Journal of Machine Learning Research}, \strong{52}, 438-448.

  Cano A and Gomez-Olmedo M and Masegosa AR and Moral S (2013). "Locally
      Averaged Bayesian Dirichlet Metrics for Learning the Structure and the
      Parameters of Bayesian Networks". \emph{International Journal of
      Approximate Reasoning}, \strong{54}, 526-540.

}
\value{

  For \code{score} with \code{by.node = TRUE}, a vector of numeric values, the
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
