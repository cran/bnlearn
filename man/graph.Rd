\name{graph utilities}
\alias{graph utilities}
\alias{acyclic}
\alias{directed}
\alias{path}
\alias{skeleton}
\alias{pdag2dag}
\alias{subgraph}
\title{ Utilities to manipulate graphs }
\description{

  Check and manipulate graph-related properties of an object of class \code{bn}.

}
\usage{
# check whether the graph is acyclic/completely directed.
acyclic(x, debug = FALSE)
directed(x)
# check whether there is a path between two nodes.
path(x, from, to, direct = TRUE, underlying.graph = FALSE,
  debug = FALSE)
# build the skeleton or a complete orientation of the graph.
skeleton(x)
pdag2dag(x, ordering)
# build a subgraph spanning a subset of nodes.
subgraph(x, nodes)
}
\arguments{
  \item{x}{an object of class \code{bn}. \code{acyclic}, \code{directed} and
       \code{path} also accept objects of class \code{bn.fit}.}
  \item{from}{a character string, the label of a node.}
  \item{to}{a character string, the label of a node (different from \code{from}).}
  \item{direct}{a boolean value. If \code{FALSE} ignore any arc between \code{from}
       and \code{to} when looking for a path.}
  \item{underlying.graph}{a boolean value. If \code{TRUE} the underlying undirected
       graph is used instead of the (directed) one from the \code{x} parameter.}
  \item{ordering}{the labels of all the nodes in the graph; their order is the
       node ordering used to set the direction of undirected arcs.}
  \item{nodes}{the labels of the nodes that induce the subgraph.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
       printed; otherwise the function is completely silent.}
}
\value{

  \code{acyclic}, \code{path} and \code{directed} return a boolean value. \cr
  \code{skeleton}, \code{pdag2dag} and \code{subgraph} return an object of
  class \code{bn}.

}
\references{

  Bang-Jensen J, Gutin G (2009). \emph{Digraphs: Theory, Algorithms and
      Applications}. Springer, 2nd edition.

}
\examples{
data(learning.test)
res = gs(learning.test)

acyclic(res)
# [1] TRUE
directed(res)
# [1] FALSE
res = pdag2dag(res, ordering = LETTERS[1:6])
res
#
#   Bayesian network learned via Constraint-based methods
#
#   model:
#     [A][C][F][B|A][D|A:C][E|B:F]
#   nodes:                                 6
#   arcs:                                  5
#     undirected arcs:                     0
#     directed arcs:                       5
#   average markov blanket size:           2.33
#   average neighbourhood size:            1.67
#   average branching factor:              0.83
#
#   learning algorithm:                    Grow-Shrink
#   conditional independence test:         Mutual Information (discrete)
#   alpha threshold:                       0.05
#   tests used in the learning procedure:  43
#   optimized:                             TRUE
#
directed(res)
# [1] TRUE
skeleton(res)
#
#   Bayesian network learned via Constraint-based methods
#
#   model:
#     [partially directed graph]
#   nodes:                                 6
#   arcs:                                  5
#     undirected arcs:                     5
#     directed arcs:                       0
#   average markov blanket size:           1.67
#   average neighbourhood size:            1.67
#   average branching factor:              0.00
#
#   learning algorithm:                    Grow-Shrink
#   conditional independence test:         Mutual Information (discrete)
#   alpha threshold:                       0.05
#   tests used in the learning procedure:  43
#   optimized:                             TRUE
#
}
\author{ Marco Scutari }
\keyword{utilities}
\keyword{graphs}
