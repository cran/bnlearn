\name{single-node local discovery}
\alias{single-node local discovery}
\alias{learn.mb}
\alias{learn.nbr}
\title{ Discover the structure around a single node }
\description{

  Learn the Markov blanket or the neighbourhood centered on a
  node.

}
\usage{
learn.mb(x, node, method, whitelist = NULL, blacklist = NULL, start = NULL,
  test = NULL, alpha = 0.05, B = NULL, debug = FALSE, optimized = TRUE)
learn.mb(x, node, method, whitelist = NULL, blacklist = NULL, start = NULL,
  test = NULL, alpha = 0.05, B = NULL, debug = FALSE, optimized = TRUE)
}
\arguments{
  \item{x}{a data frame containing the variables in the model.}
  \item{node}{a character string, the label of the node whose local 
      structure is being learned.}
  \item{method}{a character string, the label of a structure learning
     algorithm. Possible choices are \link{constraint-based algorithms}
     for \code{learn.mb} and \link{local discovery algorithms} for
     \code{learn.nbr}.}
  \item{whitelist}{a vector of character strings, the labels of the
      whitelisted nodes.}
  \item{blacklist}{a vector of character strings, the labels of the
      blacklisted nodes.}
  \item{start}{a vector of character strings, the labels of the nodes to
      be included in the Markov blanket before the learning process (in
      \code{learn.mb}). Note that the nodes in \code{start} can be removed
      from the Markov blanket by the learning algorithm, unlike the nodes
      included due to whitelisting.}
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
}
\value{

  A vector of character strings, the labels of the nodes in the Markov
  blanket (for \code{learn.mb}) or in the neighbourhood (for \code{learn.nbr}).

}
\author{ Marco Scutari }
\seealso{\link{constraint-based algorithms}, \link{local discovery algorithms}. }
\keyword{multivariate}
\keyword{models}
\keyword{graphs}
