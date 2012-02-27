\name{cpdag}
\alias{cpdag}
\alias{cextend}
\alias{vstructs}
\alias{moral}
\title{ Equivalence classes, moral graphs and consistent extenions}
\description{

  Find the equivalence class and the v-structures of a Bayesian network,
  construct its moral graph, or create a consistent extension of an
  equivalent class.

}
\usage{
cpdag(x, debug = FALSE)
cextend(x, debug = FALSE)
vstructs(x, arcs = FALSE, debug = FALSE)
moral(x, debug = FALSE)
}
\arguments{
  \item{x}{an object of class \code{bn}.}
  \item{arcs}{a boolean value. If \code{TRUE} the arcs that are part
    of at least one v-structure are returned instead of the v-structures
    themselves.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output
       is printed; otherwise the function is completely silent.}
}
\value{

  \code{cpdag} returns an object of class \code{bn}, representing the
  equivalence class. \code{moral} on the other hand returns the moral
  graph. See \code{\link{bn-class}} for details.

  \code{cextend} returns an object of class \code{bn}, representing a DAG
  that is the consistent extension of \code{x}.

  \code{vstructs} returns a matrix with either 2 or 3 columns, according
  to the value of the \code{arcs} parameter.

}
\references{

  Pearl J (2000). \emph{Causality: Models, Reasoning and Inference}.
      Cambridge University Press.

  Dor D (1992). \emph{A Simple Algorithm to Construct a Consistent
      Extension of a Partially Oriented Graph}. UCLA, Cognitive Systems
      Laboratory. Available as Technical Report R-185. 

}
\examples{
data(learning.test)
res = gs(learning.test)
cpdag(res)
#
#   Bayesian network learned via Constraint-based methods
#
#   model:
#     [partially directed graph]
#   nodes:                                 6
#   arcs:                                  5
#     undirected arcs:                     1
#     directed arcs:                       4
#   average markov blanket size:           2.33
#   average neighbourhood size:            1.67
#   average branching factor:              0.67
#
#   learning algorithm:                    Grow-Shrink
#   conditional independence test:         Mutual Information (discrete)
#   alpha threshold:                       0.05
#   tests used in the learning procedure:  43
#   optimized:                             TRUE
#
vstructs(res)
#      X   Z   Y
# [1,] "A" "D" "C"
# [2,] "B" "E" "F"
}
\author{ Marco Scutari }
\keyword{graphs}
\keyword{multivariate}
