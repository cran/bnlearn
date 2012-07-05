\name{constraint-based algorithms}
\alias{constraint-based algorithms}
\alias{gs}
\alias{iamb}
\alias{fast.iamb}
\alias{inter.iamb}
\title{ Constraint-based structure learning algorithms }
\description{

  Learn the equivalence class of a directed acyclic graph (DAG)
  from data using the Grow-Shrink (GS), the Incremental Association
  (IAMB), the Fast Incremental Association (Fast IAMB) or the
  Interleaved Incremental Association (Inter IAMB) constraint-based
  algorithms. Also use the same algorithms to learn the Markov
  blanket of a single variable.

}
\usage{
gs(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
  test = NULL, alpha = 0.05, B = NULL, debug = FALSE,
  optimized = TRUE, strict = FALSE, undirected = FALSE)
iamb(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
  test = NULL, alpha = 0.05, B = NULL, debug = FALSE,
  optimized = TRUE, strict = FALSE, undirected = FALSE)
fast.iamb(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
  test = NULL, alpha = 0.05, B = NULL, debug = FALSE,
  optimized = TRUE, strict = FALSE, undirected = FALSE)
inter.iamb(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
  test = NULL, alpha = 0.05, B = NULL, debug = FALSE,
  optimized = TRUE, strict = FALSE, undirected = FALSE)
}
\arguments{
  \item{x}{a data frame containing the variables in the model.}
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
     the learning process generate an error; otherwise they result
     in a warning.}
  \item{undirected}{a boolean value. If \code{TRUE} no attempt will be made
     to determine the orientation of the arcs; the returned (undirected)
     graph will represent the underlying structure of the Bayesian network.}
}
\value{

  An object of class \code{bn}.
  See \code{\link{bn-class}} for details.

}
\references{

  \bold{for Grow-Shrink (GS):}

  Margaritis D (2003). \emph{Learning Bayesian Network Model Structure
      from Data}. Ph.D. thesis, School of Computer Science, Carnegie-Mellon
      University, Pittsburgh, PA. Available as Technical Report CMU-CS-03-153.

  \bold{for Incremental Association (IAMB):}

  Tsamardinos I, Aliferis CF, Statnikov A (2003). "Algorithms for Large
      Scale Markov Blanket Discovery". In "Proceedings of the Sixteenth
      International Florida Artificial Intelligence Research Society
      Conference", pp. 376-381. AAAI Press.

  \bold{for Fast IAMB and Inter IAMB:}

  Yaramakala S, Margaritis D (2005). "Speculative Markov Blanket Discovery
      for Optimal Feature Selection". In "ICDM '05: Proceedings of the Fifth
      IEEE International Conference on Data Mining", pp. 809-812. IEEE
      Computer Society.

}
\author{ Marco Scutari }
\seealso{\link{local discovery algorithms},
  \link{score-based algorithms}, \link{hybrid algorithms}. }
\keyword{multivariate}
\keyword{models}
\keyword{graphs}
