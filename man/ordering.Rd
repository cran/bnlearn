\name{node ordering utilities}
\alias{node ordering utilities}
\alias{node.ordering}
\title{Partial node orderings}
\description{

  Find the partial node ordering implied by a network.

}
\usage{
node.ordering(x, debug = FALSE)
}
\arguments{
  \item{x}{an object of class \code{bn} or \code{bn.fit}.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
}
\value{

  \code{node.ordering()} returns a vector of character strings, an ordered set
  of node labels.

}
\note{

  \code{node.ordering()} supports only completely directed Bayesian networks.

}
\examples{
dag = random.graph(LETTERS[1:10])
ord = node.ordering(dag)
ord
}
\author{Marco Scutari}
\keyword{convenience functions}
\keyword{graphs}
