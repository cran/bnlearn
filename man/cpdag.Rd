\name{cpdag}
\alias{cpdag}
\alias{cextend}
\alias{vstructs}
\alias{moral}
\title{Equivalence classes, moral graphs and consistent extensions}
\description{

  Find the equivalence class and the v-structures of a Bayesian network,
  construct its moral graph, or create a consistent extension of an equivalent
  class.

}
\usage{
cpdag(x, moral = TRUE, wlbl = FALSE, debug = FALSE)
cextend(x, strict = TRUE, debug = FALSE)
vstructs(x, arcs = FALSE, moral = TRUE, debug = FALSE)
moral(x, debug = FALSE)
}
\arguments{
  \item{x}{an object of class \code{bn}.}
  \item{arcs}{a boolean value. If \code{TRUE} the arcs that are part of at least
    one v-structure are returned instead of the v-structures themselves.}
  \item{moral}{a boolean value. If \code{TRUE} we define a v-structure as in
    Pearl (2000); if \code{FALSE}, as in Koller and Friedman (2009). See below.}
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

  What kind of arc configuration is called a v-structure is not uniquely
  defined in literature. The original definition from Pearl (2000), which is
  still followed by most texts and papers, states that the two parents in the
  v-structure must not be connected by an arc. However, Koller and Friedman
  (2009) call that a \emph{immoral v-structure} and call a \emph{moral 
  v-structure} a v-structure in which the parents are linked by an arc. This
  mirrors the \emph{unshielded} versus \emph{shielded collider} naming 
  convention, but it is confusing.

  Setting \code{moral} to \code{FALSE} in \code{cpdag} and \code{vstructs} makes
  those functions follow the definition from Koller and Friedman (2009); the
  default value of \code{TRUE}, on the other hand, makes those functions follow
  the definition from Pearl (2000). The former call \emph{v-structures} both
  shielded and unshielded colliders (respectively \emph{moral v-structures} and
  \emph{immoral v-structures}); the latter requires v-structures to be 
  unshielded colliders.

  Note that arcs whose directions are dictated by the parametric assumptions of
  conditional linear Gaussian networks are preserved as directed arcs in
  \code{cpdag}.

}
\value{

  \code{cpdag} returns an object of class \code{bn}, representing the
  equivalence class. \code{moral} on the other hand returns the moral graph.
  See \code{\link{bn-class}} for details.

  \code{cextend} returns an object of class \code{bn}, representing a DAG that
  is the consistent extension of \code{x}.

  \code{vstructs} returns a matrix with either 2 or 3 columns, according to the
  value of the \code{arcs} parameter.

}
\references{

  Dor D (1992). \emph{A Simple Algorithm to Construct a Consistent Extension of
      a Partially Oriented Graph}. UCLA, Cognitive Systems Laboratory. Available
      as Technical Report R-185.

  Koller D, Friedman N (2009). \emph{Probabilistic Graphical Models: Principles
      and Techniques}. MIT Press.

  Pearl J (2009). \emph{Causality: Models, Reasoning and Inference}. Cambridge
      University Press, 2nd edition.

}
\examples{
data(learning.test)
res = gs(learning.test)
cpdag(res)
vstructs(res)
}
\author{ Marco Scutari }
\keyword{graphs}
\keyword{multivariate}
