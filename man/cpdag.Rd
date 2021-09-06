\name{cpdag}
\alias{cpdag}
\alias{cextend}
\alias{moral}
\alias{colliders}
\alias{shielded.colliders}
\alias{unshielded.colliders}
\alias{vstructs}
\title{Equivalence classes, moral graphs and consistent extensions}
\description{

  Find the equivalence class and the v-structures of a Bayesian network,
  construct its moral graph, or create a consistent extension of an equivalent
  class.

}
\usage{
cpdag(x, wlbl = FALSE, debug = FALSE)
cextend(x, strict = TRUE, debug = FALSE)
moral(x, debug = FALSE)

colliders(x, arcs = FALSE, debug = FALSE)
shielded.colliders(x, arcs = FALSE, debug = FALSE)
unshielded.colliders(x, arcs = FALSE, debug = FALSE)
vstructs(x, arcs = FALSE, debug = FALSE)
}
\arguments{
  \item{x}{an object of class \code{bn} or \code{bn.fit} (with the exception of
    \code{cextend}, which only accepts objects of class \code{bn}).}
  \item{arcs}{a boolean value. If \code{TRUE} the arcs that are part of at least
    one v-structure are returned instead of the v-structures themselves.}
  \item{wlbl}{a boolean value. If \code{TRUE} arcs whose directions have been
    fixed by a whitelist or a by blacklist are preserved when constructing
    the CPDAG.}
  \item{strict}{a boolean value. If no consistent extension is possible and
    \code{strict} is \code{TRUE}, an error is generated; otherwise a partially
    extended graph is returned with a warning.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
}
\details{

  Note that arcs whose directions are dictated by the parametric assumptions of
  conditional linear Gaussian networks are preserved as directed arcs in
  \code{cpdag()}.

}
\value{

  \code{cpdag()} returns an object of class \code{bn}, representing the
  equivalence class. \code{moral} on the other hand returns the moral graph.
  See \code{\link{bn-class}} for details.

  \code{cextend()} returns an object of class \code{bn}, representing a DAG that
  is the consistent extension of \code{x}.

  \code{vstructs()} returns a matrix with either 2 or 3 columns, according to the
  value of the \code{arcs} argument.

}
\references{

  Dor D (1992). \emph{A Simple Algorithm to Construct a Consistent Extension of
    a Partially Oriented Graph}. UCLA, Cognitive Systems Laboratory.

  Koller D, Friedman N (2009). \emph{Probabilistic Graphical Models: Principles
    and Techniques}. MIT Press.

  Pearl J (2009). \emph{Causality: Models, Reasoning and Inference}. Cambridge
    University Press, 2nd edition.

}
\examples{
data(learning.test)
res = hc(learning.test)
cpdag(res)
vstructs(res)
}
\author{Marco Scutari}
\keyword{graphs}
\keyword{structure learning}
