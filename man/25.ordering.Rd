\name{node ordering tools}
\alias{node.ordering}
\alias{ordering2blacklist}
\title{ Partial node ordering utilities }
\description{

  Detect the partial node ordering implied by a network or
  generate the blacklist implied by a complete node ordering.

}
\usage{
  node.ordering(x, debug = FALSE)

  ordering2blacklist(nodes)
}
\arguments{
  \item{x}{an object of class \code{bn}.}
  \item{nodes}{a vector of character strings, the labels of the nodes.
    The ordering of the labels must reflect the partial node ordering
    of the nodes in the graph.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
       printed; otherwise the function is completely silent.}
}
\value{

  \code{node.ordering} return a vector of character strings, an ordered
  set of node labels.

  \code{ordering2blacklist} returns a sanitized \code{blacklist} (a
  two-column matrix, whose columns are labeled \code{from} and \code{to}).

}
\note{

  \code{node.ordering} and \code{ordering2blacklist} support only
  completely directed Bayesian networks.

}
\examples{
data(learning.test)
res = gs(learning.test, optimized = TRUE)
res$learning$ntests
# [1] 43
res = set.arc(res, "A", "B")
ord = node.ordering(res)
ord
# [1] "A" "C" "F" "B" "D" "E"

## partial node ordering saves us two tests in the v-structure
## detection step of the algorithm.
gs(learning.test, blacklist = ordering2blacklist(ord))$learning$ntests
# [1] 41
}
\author{ Marco Scutari }
\keyword{utilities}
\keyword{multivariate}
