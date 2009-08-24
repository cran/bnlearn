\name{mmpc}
\alias{mmpc}
\title{ Max-Min Parents and Children (MMPC) learning algorithm }
\description{

  Estimate the underlying structure of a directed acyclic graph (DAG)
  from data using the Max-Min Parents and Children (MMPC)
  constraint-based algorithm.

}
\usage{
  mmpc(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
    test = NULL, alpha = 0.05, B = NULL, debug = FALSE,
    optimized = TRUE, strict = FALSE)
}
\arguments{
  \item{x}{a data frame, containing the variables in the model.}
  \item{cluster}{an optional cluster object from package \pkg{snow}.
      See \code{\link{snow integration}} for details and a simple
      example.}
  \item{whitelist}{a data frame with two columns (optionally labeled
     "from" and "to"), containing a set of arcs to be included in the
      graph. }
  \item{blacklist}{a data frame with two columns (optionally labeled
     "from" and "to"), containing a set of arcs not to be included in
      the graph. }
  \item{test}{a character string, the label of the conditional
      independence test to be used in the algorithm. If none is
      specified, the default test statistic is the \emph{mutual information}
      for discrete data sets and the \emph{linear correlation} for
      continuous ones. See \code{\link{bnlearn-package}} for details.}
  \item{alpha}{a numeric value, the target nominal type I error rate.}
  \item{B}{a positive integer, the number of permutations considered
      for each permutation test. It will be ignored with a warning if
      the conditional independence test specified by the \code{test}
      argument is not a permutation test.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output
      is printed; otherwise the function is completely silent.}
  \item{optimized}{a boolean value. See \code{\link{bnlearn-package}}
      for details.}
  \item{strict}{a boolean value. If \code{TRUE} conflicting results in
      the learning process generate an error; otherwise they result in
      a warning.}
}
\value{

  An object of class \code{bn}.
  See \code{\link{bn-class}} for details.

}
\references{

  Tsamardinos I, Aliferis CF, Statnikov A (2003). "Time and Sample
     Efficient Discovery of Markov Blankets and Direct Causal Relations".
     In "KDD '03: Proceedings of the Ninth ACM SIGKDD International
     Conference on Knowledge Discovery and Data Mining", pp. 673-678. ACM.

  Tsamardinos I, Brown LE, Aliferis CF (2006). "The Max-Min Hill-Climbing
      Bayesian Network Structure Learning Algorithm". \emph{Machine Learning},
      \strong{65}(1), 31-78.

}
\author{ Marco Scutari }
\seealso{\code{\link{constraint-based algorithms}}, \code{\link{hc}},
  \code{\link{hybrid algorithms}}. }
\keyword{multivariate}
\keyword{models}
\keyword{graphs}
