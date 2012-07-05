\name{dsep}
\alias{dsep}
\title{ Test d-separation }
\description{

  Check whether two nodes are d-separated.

}
\usage{
dsep(bn, x, y, z)
}
\arguments{
  \item{bn}{an object of class \code{bn}.}
  \item{x,y}{a character string, the label of a node.}
  \item{z}{an optional vector of character strings, the label of the
    (candidate) d-separating nodes. It defaults to the empty set.}
}
\value{

  \code{dsep} returns \code{TRUE} if \code{x} and \code{y} are
  d-separated by \code{z}, and \code{FALSE} otherwise.

}
\references{

  Koller D, Friedman N (2009). \emph{Probabilistic Graphical Models:
      Principles and Techniques}. MIT Press.

}
\examples{
bn = model2network("[A][C|A][B|C]")
dsep(bn, "A", "B", "C")
# [1] TRUE
bn = model2network("[A][C][B|A:C]")
dsep(bn, "A", "B", "C")
# [1] FALSE
}
\author{ Marco Scutari }
\keyword{graphs}
\keyword{multivariate}
