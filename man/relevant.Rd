\name{relevant}
\alias{relevant}
\title{Identify Relevant Nodes Without Learning the Bayesian network}
\description{

  Identify all the nodes relevant to compute all the conditional probability
  distributions for a given set of nodes.

}
\usage{
  relevant(target, context, data, test, alpha, B, debug = FALSE)
}
\arguments{
  \item{target}{a vector of character strings, the labels of nodes whose
    conditional probability distributions are of interest.}
  \item{context}{a vector of character strings, the labels of nodes on
    which to condition the independence tests.}
  \item{data}{a data frame containing either numeric or factor columns.}
  \item{test}{a character string, the label of the conditional independence
    test to be used in the algorithm. If none is specified, the default test
    statistic is the \emph{mutual information} for categorical variables, the
    Jonckheere-Terpstra test for ordered factors and the \emph{linear
    correlation} for continuous variables. See \code{\link{bnlearn-package}}
    for details.}
  \item{alpha}{a numeric value, the target nominal type I error rate. If none
    is specified, the default value is \code{0.05}.}
  \item{B}{a positive integer, the number of permutations considered for each
    permutation test. It will be ignored with a warning if the conditional
    independence test specified by the \code{test} argument is not a
    permutation test.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
}
\note{

  This algorithms selects all the nodes that are relevant at all, not only
  those that are significantly so. Therefore, to be discarded a node must
  be completely unrelated to any of the \code{target} nodes, not just weakly
  dependent. On the good side, relevant nodes are correctly identified even
  for data sets whose probability structure is not faithful to any directed
  acyclic graph.

}
\value{

  \code{relevant} returns a vector of character strings, the labels of the
  relevant nodes.

}
\examples{
data(learning.test)
X = as.factor(sample(c("x1", "x2"), nrow(learning.test), replace = TRUE))
relevant("A", data = cbind(learning.test, X))
relevant("A", context = "B", data = learning.test,)
}
\references{

  Pena JM, Nilsson R, Bjorkegren J, Tegner J (2006). "Identifying the Relevant
      Nodes Without Learning the Model". In "Proceedings of the 22nd Conference
      on Uncertainty in Artificial Intelligence (UAI2006)", pp. 367-374.

}
\author{ Marco Scutari }
\keyword{multivariate}
\keyword{models}
