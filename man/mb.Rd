\name{misc utilities}
\alias{misc utilities}
\alias{mb}
\alias{nbr}
\alias{arcs}
\alias{arcs<-}
\alias{directed.arcs}
\alias{undirected.arcs}
\alias{narcs}
\alias{nodes}
\alias{amat}
\alias{amat<-}
\alias{parents}
\alias{parents<-}
\alias{children}
\alias{children<-}
\alias{root.nodes}
\alias{leaf.nodes}
\alias{nparams}
\title{ Miscellaneous utilities }
\description{

  Assign or extract various quantities of interest from an object of
  class \code{bn} of \code{bn.fit}.

}
\usage{
## nodes
nodes(x)
mb(x, node)
nbr(x, node)
parents(x, node)
parents(x, node, debug = FALSE) <- value
children(x, node)
children(x, node, debug = FALSE) <- value
root.nodes(x)
leaf.nodes(x)

## arcs
arcs(x)
arcs(x, ignore.cycles = FALSE, debug = FALSE) <- value
directed.arcs(x)
undirected.arcs(x)
narcs(x)

## adjacency matrix
amat(x)
amat(x, ignore.cycles = FALSE, debug = FALSE) <- value

## graphs
nparams(x, data, debug = FALSE)
}
\arguments{
  \item{x}{an object of class \code{bn} or \code{bn.fit}. The replacement
       form of \code{parents}, \code{children}, \code{arcs} and \code{amat}
       require an object of class \code{bn}.}
  \item{node}{a character string, the label of a node.}
  \item{value}{either a vector of character strings (for \code{parents} and
       \code{children}), an adjacency matrix (for \code{amat}) or a data
       frame with two columns (optionally labeled "from" and "to", for
       \code{arcs}).}
  \item{data}{a data frame containing the data the Bayesian network
       was learned from. It's only needed if \code{x} is an object of
       class \code{bn}.}
  \item{ignore.cycles}{a boolean value. If \code{TRUE} the returned network
       will not be checked for cycles.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
       printed; otherwise the function is completely silent.}
}
\details{

  The number of parameters of a discrete Bayesian network is defined as the
  sum of the number of logically independent parameters of each node given its
  parents (Chickering, 1995). For Gaussian Bayesian networks the distribution
  of each node can be viewed as a linear regression, so it has a number of
  parameters equal to the number of the parents of the node plus one (the
  intercept) as per Neapolitan (2003).

}
\value{

  \code{mb}, \code{nbr}, \code{nodes}, \code{parents}, \code{root.nodes}
  and \code{leaf.nodes} return a vector of character strings. \cr
  \code{arcs} returns a matrix of two columns of character strings. \cr
  \code{narcs} returns the number of arcs in the graph.
  \code{amat} returns a matrix of 0/1 integer values. \cr
  \code{nparams} returns an integer. \cr

}
\references{

  Chickering DM (1995). "A Transformational Characterization of Equivalent
      Bayesian Network Structures". In "UAI '95: Proceedings of the Eleventh
      Annual Conference on Uncertainty in Artificial Intelligence", pp. 87-98.
      Morgan Kaufmann.

  Neapolitan RE (2003). \emph{Learning Bayesian Networks}. Prentice Hall.

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
root.nodes(res)
# [1] "C" "F"
## the leaf nodes of the graph.
leaf.nodes(res)
# [1] "D" "E"
## number of parameters of the Bayesian network.
res = set.arc(res, "A", "B")
nparams(res, learning.test)
# [1] 41
}
\author{ Marco Scutari }
\keyword{utilities}
\keyword{graphs}
