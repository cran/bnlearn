\name{node ordering utilities}
\alias{node ordering utilities}
\alias{node.ordering}
\alias{ordering2blacklist}
\alias{tiers2blacklist}
\title{ Utilities dealing with partial node orderings }
\description{

  Find the partial node ordering implied by a network or
  generate the blacklist implied by a complete node ordering.

}
\usage{
node.ordering(x, debug = FALSE)
ordering2blacklist(nodes)
tiers2blacklist(nodes)
}
\arguments{
  \item{x}{an object of class \code{bn} or \code{bn.fit}.}
  \item{nodes}{a node ordering, see below.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
       printed; otherwise the function is completely silent.}
}
\details{

  \code{ordering2blacklist} takes a vector of character strings (the
  labels of the nodes), which specifies a complete node ordering. An
  object of class \code{bn} or \code{bn.fit}; in that case, the node
  ordering is derived by the graph. In both cases, the blacklist returned
  by \code{ordering2blacklist} contains all the possible arcs that
  violate the specified node ordering.

  \code{tiers2blacklist} takes (again) a vector of character strings
  (the labels of the nodes), which specifies a complete node ordering,
  or a list of character vectors, which specifies a partial node
  ordering. In the latter case, all arcs going from a node in a
  particular element of the list (sometimes known as \emph{tier}) to
  a node in one of the previous elements are blacklisted. Arcs between
  nodes in the same element are not blacklisted.

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

tiers2blacklist(list(LETTERS[1:3], LETTERS[4:6]))
#       from to 
#  [1,] "D"  "A"
#  [2,] "E"  "A"
#  [3,] "F"  "A"
#  [4,] "D"  "B"
#  [5,] "E"  "B"
#  [6,] "F"  "B"
#  [7,] "D"  "C"
#  [8,] "E"  "C"
#  [9,] "F"  "C"
}
\author{ Marco Scutari }
\keyword{utilities}
\keyword{multivariate}
