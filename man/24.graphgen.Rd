\name{graph generation tools}
\alias{empty.graph}
\alias{random.graph}
\title{ Generate an empty or random graph }
\description{
  Generate an empty or random graph from a given set of nodes.
}
\usage{
  empty.graph(nodes)
  random.graph(nodes, prob = 0.5)
}
\arguments{
  \item{nodes}{a vector of character strings, the labels of the nodes.}
  \item{prob}{a numeric value, the probability of connecting a node to
    another node with higher ordering.}
}
\value{

  Both \code{empty.graph} and \code{random.graph} return an
  object of class \code{bn}.

}
\note{

  The \code{random.graph} function uses the same algorithm as the
  \code{randomDAG} function in package \pkg{pcalg} and generates
  graphs whose node ordering is given by the order of the labels
  in the \code{nodes} parameter.

}
\examples{
empty.graph(LETTERS[1:8])
#
#  Bayesian network learned via Conditional Independence methods
#
#  model:
#    [A][B][C][D][E][F][G][H]
#  nodes:                                 8
#  arcs:                                  0
#    undirected arcs:                     0
#    directed arcs:                       0
#  average markov blanket size:           0.00
#  average neighbourhood size:            0.00
#  average branching factor:              0.00
#
#  learning algorithm:                    random/generated
#  tests used in the learning procedure:  0
#
random.graph(LETTERS[1:8])
# <insert the description of a random graph here>
plot(random.graph(LETTERS[1:8]))
plot(random.graph(LETTERS[1:8], prob = 0.2))
}
\author{ Marco Scutari }
\keyword{utilities}
\keyword{graphs}
