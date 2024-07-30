\name{single-node local discovery}
\alias{single-node local discovery}
\alias{learn.mb}
\alias{learn.nbr}
\title{Discover the structure around a single node}
\description{

  Learn the Markov blanket or the neighbourhood centered on a node.

}
\usage{
learn.mb(x, node, method, whitelist = NULL, blacklist = NULL, start = NULL,
  test = NULL, alpha = 0.05, ..., max.sx = NULL, debug = FALSE)
learn.nbr(x, node, method, whitelist = NULL, blacklist = NULL,
  test = NULL, alpha = 0.05, ..., max.sx = NULL, debug = FALSE)
}
\arguments{
  \item{x}{a data frame containing the variables in the model.}
  \item{node}{a character string, the label of the node whose local structure
    is being learned.}
  \item{method}{a character string, the label of a structure learning algorithm.
    Possible choices are listed in \link{structure learning}.}
  \item{whitelist}{a vector of character strings, the labels of the whitelisted
    nodes.}
  \item{blacklist}{a vector of character strings, the labels of the blacklisted
    nodes.}
  \item{start}{a vector of character strings, the labels of the nodes to be
    included in the Markov blanket before the learning process (in
    \code{learn.mb}). Note that the nodes in \code{start} can be removed from
    the Markov blanket by the learning algorithm, unlike the nodes included due
    to whitelisting.}
  \item{test}{a character string, the label of the conditional independence test
    to be used in the algorithm. If none is specified, the default test
    statistic is the \emph{mutual information} for categorical variables, the
    Jonckheere-Terpstra test for ordered factors and the \emph{linear
    correlation} for continuous variables. See \code{\link{independence tests}}
    for details.}
  \item{alpha}{a numeric value, the target nominal type I error rate.}
  \item{\dots}{optional arguments to be passed to the test specified by
    \code{test}. See \code{\link{ci.test}} for details.}
  \item{max.sx}{a positive integer, the maximum allowed size of the conditioning
    sets used in conditional independence tests. The default is that there is
    no limit on size.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
}
\note{

  All algorithms used by \code{learn.mb()} and \code{learn.nbr()} accept
  incomplete data, which they handle by computing individual conditional
  independence tests on locally complete observations.

}
\value{

  A vector of character strings, the labels of the nodes in the Markov blanket
  (for \code{learn.mb()}) or in the neighbourhood (for \code{learn.nbr()}).

}
\examples{
learn.mb(learning.test, node = "D", method = "iamb")
learn.mb(learning.test, node = "D", method = "iamb", blacklist = c("A", "F"))

learn.nbr(gaussian.test, node = "F", method = "si.hiton.pc", whitelist = "D")
}
\author{Marco Scutari}
\seealso{\link{constraint-based algorithms}.}
\keyword{local learning}
