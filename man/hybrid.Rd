\name{hybrid algorithms}
\alias{hybrid algorithms}
\alias{rsmax2}
\alias{mmhc}
\title{ Hybrid structure learning algorithms }
\description{

  Learn the structure of a Bayesian network with the Max-Min
  Hill Climbing (MMHC) and the more general 2-phase Restricted
  Maximization (RSMAX2) hybrid algorithms.

}
\usage{
rsmax2(x, whitelist = NULL, blacklist = NULL, restrict,
  maximize = "hc", test = NULL, score = NULL, alpha = 0.05,
  B = NULL, ..., maximize.args = list(), optimized = TRUE,
  strict = FALSE, debug = FALSE)
mmhc(x, whitelist = NULL, blacklist = NULL, test = NULL,
  score = NULL, alpha = 0.05, B = NULL, ..., restart = 0,
  perturb = 1, max.iter = Inf, optimized = TRUE,
  strict = FALSE, debug = FALSE)
}
\arguments{
  \item{x}{a data frame containing the variables in the model.}
  \item{whitelist}{a data frame with two columns (optionally labeled
     "from" and "to"), containing a set of arcs to be included in the
     graph. }
  \item{blacklist}{a data frame with two columns (optionally labeled
     "from" and "to"), containing a set of arcs not to be included in
     the graph. }
  \item{restrict}{a character string, the constraint-based algorithm
     to be used in the \dQuote{restrict} phase. Possible values are
     \code{gs}, \code{iamb}, \code{fast.iamb}, \code{inter.iamb} and
     \code{mmpc}. See \code{\link{bnlearn-package}} and the
     documentation of each algorithm for details.}
  \item{maximize}{a character string, the score-based algorithm
     to be used in the \dQuote{maximize} phase. Possible values are
     \code{hc} and \code{tabu}. See \code{\link{bnlearn-package}} for
     details.}
  \item{test}{a character string, the label of the conditional
     independence test to be used by the constraint-based algorithm.
     If none is specified, the default test statistic is the
     \emph{mutual information} for discrete data sets and the
     \emph{linear correlation} for continuous ones. See
     \code{\link{bnlearn-package}} for details.}
  \item{score}{a character string, the label of the network score to
      be used in the score-based algorithm. If none is specified, the
      default score is the \emph{Bayesian Information Criterion} for
      both discrete and continuous data sets. See \code{\link{bnlearn-package}}
      for details.}
  \item{alpha}{a numeric value, the target nominal type I error rate of
     the conditional independence test.}
  \item{B}{a positive integer, the number of permutations considered
     for each permutation test. It will be ignored with a warning if
     the conditional independence test specified by the \code{test}
     argument is not a permutation test.}
  \item{\dots}{additional tuning parameters for the network score used
     by the score-based algorithm. See \code{\link{score}} for details.}
  \item{maximize.args}{a list of arguments to be passed to the score-based
     algorithm specified by \code{maximize}, such as \code{restart} for
     hill-climbing or \code{tabu} for tabu search.}
   \item{restart}{an integer, the number of random restarts for the
      score-based algorithm.}
   \item{perturb}{an integer, the number of attempts to randomly
      insert/remove/reverse an arc on every random restart.}
   \item{max.iter}{an integer, the maximum number of iterations for the
      score-based algorithm.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output
     is printed; otherwise the function is completely silent.}
  \item{optimized}{a boolean value. See \code{\link{bnlearn-package}}
     for details.}
  \item{strict}{a boolean value. If \code{TRUE} conflicting results in
     the learning process generate an error; otherwise they result
     in a warning.}
}
\note{

  \code{mmhc} is simply \code{rshc} with \code{restrict} set to
  \code{mmpc} and \code{maximize} set to \code{hc}.

}
\value{

  An object of class \code{bn}.
  See \code{\link{bn-class}} for details.

}
\references{

  Tsamardinos I, Brown LE, Aliferis CF (2006). "The Max-Min Hill-Climbing
      Bayesian Network Structure Learning Algorithm". \emph{Machine Learning},
      \strong{65}(1), 31-78.

}
\author{ Marco Scutari }
\seealso{\link{local discovery algorithms},
  \link{score-based algorithms}, \link{constraint-based algorithms}. }
\keyword{multivariate}
\keyword{models}
\keyword{graphs}
