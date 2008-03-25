\name{network tools}
\alias{mb}
\alias{nbr}
\alias{arcs}
\alias{arcs<-}
\alias{nodes}
\alias{amat}
\alias{amat<-}
\alias{parents}
\alias{parents<-}
\alias{children}
\alias{children<-}
\alias{rootnodes}
\alias{leafnodes}
\alias{nparams}
\alias{acyclic}
\alias{path}
\title{ Utilities to manipulate graphs }
\description{

 Assign or extract various quantities of interest from an object of class \code{bn}.

}
\usage{
  ## nodes
  nodes(x)
  mb(x, node, rebuild = FALSE)
  nbr(x, node, rebuild = FALSE)
  parents(x, node, rebuild = FALSE)
  parents(x, node, debug = FALSE) <- value
  children(x, node, rebuild = FALSE)
  children(x, node, debug = FALSE) <- value
  rootnodes(x)
  leafnodes(x)

  ## arcs
  arcs(x)
  arcs(x, debug = FALSE) <- value

  ## adjacency matrix
  amat(x)
  amat(x, debug = FALSE) <- value

  ## graphs
  nparams(x, data, debug = FALSE)
  acyclic(x, directed, debug = FALSE)
  path(x, from, to, direct = TRUE)
}
\arguments{
  \item{x}{an object of class "bn".}
  \item{node}{a character string, the label of a node.}
  \item{from}{a character string, the label of a node.}
  \item{to}{a character string, the label of a node (different from \code{from}).}
  \item{value}{either a vector of character strings (for \code{parents} and
       \code{children}), an adjacency matrix (for \code{amat}) or a data
       frame with two columns (optionally labeled "from" and "to", for
       \code{arcs})}.
  \item{data}{a data frame, containing the data the Bayesian network
       was learned from.}
  \item{rebuild}{a boolean value. If \code{TRUE} the return value is rebuilt
       from scratch using the arc list; otherwise the cached value are returned.}
  \item{direct}{a boolean value. If \code{FALSE} ignore any arc between \code{from}
       and \code{to} when looking for a path.}
  \item{directed}{a boolean value. If \code{TRUE} the graph is assumed to be
       completely directed (no undirected arcs), and a faster cycle detection
       algorithm is used.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
       printed; otherwise the function is completely silent.}
}
\details{

  The number of parameters of a discrete Bayesian network is computed as the
  sum of the number of logically independent parameters of each node given its
  parents (Chickering, 1995).

}
\value{

  \code{mb}, \code{nbr}, \code{nodes}, \code{parents}, \code{rootnodes}
  and \code{leafnodes} return a vector of character strings.
  \code{arcs} returns a matrix of two columns of character strings.
  \code{amat} returns a matrix of 0/1 numeric values.
  \code{nparams} returns an integer.
  \code{acyclic} and \code{path} return a boolean value.

}
\note{

  \code{nparams} supports only completely directed discrete Bayesian networks.

}
\references{

  D. M. Chickering. A Transformational Characterization of Equivalent Bayesian
    Network Structures. In Proceedins of 11th Conference on Uncertainty in
    Artificial Intelligence, pages 87-98. Morgan Kaufmann Publishers Inc., 1995.

}
\examples{
data(learning.test)
res = gs(learning.test)

##  the Markov blanket of A.
mb(res, "A")
# [1] "B" "D" "C"
## the neighbourhood of F.
nbr(res, "F")
# [1] "E"
## the arcs in the graph.
arcs(res)
#      from to
# [1,] "A"  "B"
# [2,] "A"  "D"
# [3,] "B"  "A"
# [4,] "B"  "E"
# [5,] "C"  "D"
# [6,] "F"  "E"
## the nodes of the graph.
nodes(res)
# [1] "A" "B" "C" "D" "E" "F"
## the adjacency matrix for the nodes of the graph.
amat(res)
#   A B C D E F
# A 0 1 0 1 0 0
# B 1 0 0 0 1 0
# C 0 0 0 1 0 0
# D 0 0 0 0 0 0
# E 0 0 0 0 0 0
# F 0 0 0 0 1 0
## the parents of D.
parents(res, "D")
# [1] "A" "C"
## the children of A.
children(res, "A")
# [1] "D"
## the root nodes of the graph.
rootnodes(res)
# [1] "C" "F"
## the leaf nodes of the graph.
leafnodes(res)
# [1] "D" "E"
## number of parameters of the Bayesian network.
res = set.arc(res, "A", "B")
nparams(res, learning.test)
# [1] 41
}
\author{ Marco Scutari }
\keyword{utilities}
\keyword{graphs}
