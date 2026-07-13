\name{graph utilities}
\alias{graph utilities}
\alias{acyclic}
\alias{directed}
\alias{valid.dag}
\alias{valid.cpdag}
\alias{valid.ug}
\alias{path.exists}
\alias{path}
\alias{path,bn-method}
\alias{path,bn.fit-method}
\alias{path,bn.naive-method}
\alias{path,bn.tan-method}
\alias{skeleton}
\alias{pdag2dag}
\alias{subgraph}
\alias{perturb}
\title{Utilities to manipulate graphs}
\description{

  Check and manipulate graph-related properties of an object of class \code{bn}.

}
\usage{
# check whether the graph is acyclic/completely directed.
acyclic(x, directed = FALSE, debug = FALSE)
directed(x)
valid.dag(x, debug = FALSE)
valid.cpdag(x, debug = FALSE)
valid.ug(x, debug = FALSE)
# check whether there is a path between two nodes.
path.exists(x, from, to, direct = TRUE, underlying.graph = FALSE, debug = FALSE)
# build the skeleton or a complete orientation of the graph.
skeleton(x)
pdag2dag(x, ordering)
# build a subgraph spanning a subset of nodes.
subgraph(x, nodes)
# perturb a network by random arc additions, removals and reversals.
perturb(x, nops, ops = c("set", "drop", "reverse"), maxp = Inf, debug = FALSE)
}
\arguments{
  \item{x}{an object of class \code{bn}. \code{skeleton()}, \code{acyclic()},
    \code{directed()} and \code{path.exists()} also accept objects of class
    \code{bn.fit}.}
  \item{from}{a character string, the label of a node.}
  \item{to}{a character string, the label of a node (different from
    \code{from}).}
  \item{direct}{a boolean value. If \code{FALSE}, ignore any arc between
    \code{from} and \code{to} when looking for a path.}
  \item{underlying.graph}{a boolean value. If \code{TRUE}, the underlying
    undirected graph is used instead of the (directed) one from the \code{x}
    argument.}
  \item{ordering}{the labels of all the nodes in the graph; their order is the
    node ordering used to set the direction of undirected arcs.}
  \item{nodes}{the labels of the nodes that induce the subgraph.}
  \item{directed}{a boolean value. If \code{TRUE}, only completely directed
    cycles are considered. Otherwise, undirected arcs will also be considered
    and treated as arcs present in both directions.}
  \item{nops}{a positive integer value, the number of arc operations applied to
    the network to perturb it.}
  \item{ops}{a vector of character strings, one or more of \code{"set"},
    \code{"drop"} or \code{"reverse"}. It defines which arc operations are
    applied to the network to perturb it.}
  \item{maxp}{the maximum number of parents allowed for a node after each arc
    operation. The default value is \code{Inf}.}
  \item{debug}{a boolean value. If \code{TRUE}, a lot of debugging output is
    printed. Otherwise, the function is completely silent.}
}
\value{

  \code{acyclic()}, \code{path()} and \code{directed()} return a boolean value. \cr
  \code{skeleton()}, \code{pdag2dag()}, \code{subgraph()} and \code{perturb()}
  return an object of class \code{bn}.

  \code{valid.dag()}, \code{valid.cpdag()} and \code{valid.ug()} return a
  boolean value.

}
\references{

  Bang-Jensen J, Gutin G (2009). \emph{Digraphs: Theory, Algorithms and
    Applications}. Springer, 2nd edition.

  Andersson SA, Madigan D, Perlman MD (1997). "A Characterization of Markov
    Equivalence Classes for Acyclic Digraphs." The Annals of Statistics,
    \strong{25}(2):505--541.

}
\examples{
data(learning.test)
cpdag = pc.stable(learning.test)

acyclic(cpdag)
directed(cpdag)
dag = pdag2dag(cpdag, ordering = LETTERS[1:6])
dag
directed(dag)
skeleton(dag)
perturb(dag, nops = 3, debug = TRUE)
}
\author{Marco Scutari}
\keyword{graphs}
