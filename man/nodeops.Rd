\name{node operations}
\alias{node operations}
\alias{add.node}
\alias{remove.node}
\alias{rename.nodes}
\alias{nodes<-}
\alias{nodes<-,bn-method}
\alias{nodes<-,bn.fit-method}
\alias{nodes<-,bn.naive-method}
\alias{nodes<-,bn.tan-method}
\title{Manipulate nodes in a graph}
\description{

  Add, remove and rename nodes in a graph.

}
\usage{
# add and remove nodes.
add.node(x, node)
remove.node(x, node)

# re-label nodes.
rename.nodes(x, names)
\S4method{nodes}{bn}(object) <- value
\S4method{nodes}{bn.fit}(object) <- value
}
\arguments{
  \item{x}{an object of class \code{bn} for \code{add.node()} and
    \code{remove.node()}; an object of class \code{bn} or \code{bn.fit} for
    \code{rename.nodes()}.}
  \item{object}{an object of class \code{bn} or \code{bn.fit}.}
  \item{node}{a character string, the label of a node.}
  \item{value, names}{a vector of character strings, the new set of labels
    that wll be used as to rename the nodes.}
}
\details{

  \code{add.node()} adds a new (isolated) node to an existing \code{bn} object.

  \code{remove.node()} removes a node from a \code{bn} object.

  \code{rename.nodes()} replaces the node labels with new ones, relabelling the
  whole node set. The assignment method for \code{nodes()} is an alias of
  \code{rename.nodes()}.

}
\value{

  \code{add.node()}, \code{remove.node()} and \code{rename.nodes()} return an
  updated \code{bn} object.

}
\examples{
dag = random.graph(LETTERS[1:5])
add.node(dag, "Z")
remove.node(dag, "A")

nodes(dag)
nodes(dag) = LETTERS[6:10]
nodes(dag)
}
\author{Marco Scutari}
\keyword{graphs}
