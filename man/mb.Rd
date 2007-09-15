\name{mb}
\alias{mb}
\alias{nbr}
\alias{arcs}
\alias{nodes}
\alias{amat}
\alias{parents}
\title{ Utility functions to manipulate graphs }
\description{

 Extract various quantities of interest from an object of class \code{bn}.

}
\usage{

  mb(x, node)
  nbr(x, node)
  arcs(x)
  nodes(x)
  amat(x)
  parents(x, node)

}
\arguments{
  \item{x}{an object of class "bn".}
  \item{node}{a character string, the label of a node.}
}
\value{

  \code{mb}, \code{nbr}, \code{nodes} and \code{parents} return an array of character strings.
  \code{arcs} returns a matrix of two columns of character strings.
  \code{amat} returns a matrix of 0/1 numeric values.

}
\examples{

data(learning.test)
res = gs(learning.test)

# the Markov blanket of A.
mb(res, "A")
# [1] "B" "D" "C"
# the neighbourhood of F.
nbr(res, "F")
# [1] "B" "E"
# the arcs in the graph.
arcs(res)
#      from to 
# [1,] "A"  "D"
# [2,] "B"  "A"
# [3,] "C"  "D"
# [4,] "E"  "B"
# [5,] "E"  "F"
# [6,] "F"  "B"
# [7,] "F"  "E"
# the nodes of the graph.
nodes(res)
# [1] "A" "B" "C" "D" "E" "F"
# the adjacency matrix for the nodes of the graph.
amat(res)
# A 0 0 0 1 0 0
# B 1 0 0 0 0 0
# C 0 0 0 1 0 0
# D 0 0 0 0 0 0
# E 0 1 0 0 0 1
# F 0 1 0 0 1 0
# the parents of D.
parents(res, "D")
# [1] "A" "C"

}
\author{ Marco Scutari }
\keyword{utilities}
