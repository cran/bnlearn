\name{ctsdag}
\alias{ctsdag}
\title{Equivalence classes in the presence of interventions}
\description{
  Generate the partially directed acyclic graph representing the equivalence
  class of a Bayesian network learned using interventions.
}
\usage{
  ctsdag(x, exp, learning = FALSE, debug = FALSE)
}
\arguments{
  \item{x}{an object of class \code{bn}, the network from which to compute
    the PDAG.}
  \item{exp}{a vector of character strings, the labels of the node that are the
    targets of the interventions. If no targets are provided, \code{ctsdag()}
    just reverts to \code{cpdag()}.}
  \item{learning}{a boolean value. If \code{TRUE}, interventions, whitelists
    and blacklists used in learning the structure of \code{x} will be taken
    into account in contructing the PDAG. These interventions will be applied
    in addition to those provided via the \code{exp} argument.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
}
\details{

  \code{ctsdag()} extends \code{cpdag()} by incorporating interventions in
  constructing the partially directed acyclic graph that represents the
  equivalence class of \code{x}; it preserves the directions of arcs that
  are compelled because they are incident on the target nodes specified by the
  \code{exp} argument. This assumes do-calculus model of targeted
  interventions with no unknown side-effects.

  It also takes into account prior arc probabilities used in structure learning,
  ensuring that DAGs are equivalent in posterior probability only if they are
  equivalent in prior probability. This is not the case for graph priors other
  than the uniform (\code{uniform}) and marginal uniform priors
  (\code{marginal}, see \code{\link{bn-class}} for details).

}
\value{
  \code{ctsdag} returns an object of class \code{bn}, representing the
  equivalence class. See \code{\link{bn-class}} for details.
}
\references{

  Castelo R, Siebes A (2000). "Priors on Network Structures. Biasing the Search
    for Bayesian Networks". \emph{International Journal of Approximate
    Reasoning}, \strong{24}(1):39--57.

  Chickering DM (1995). "A Transformational Characterization of Equivalent
	Bayesian Network Structures". \emph{Proceedings of the Eleventh Annual
	Conference on Uncertainty in Artificial Intelligence}, 87--98.

  Ness RO, Sachs K, Mallick P, Vitek O (2017). "A Bayesian Active Learning
    Experimental Design for Inferring Signaling Networks". \emph{International
    Conference on Research in Computational Molecular Biology}, 134--156.

  Tian J, Pearl J (2001). "Causal Discovery from Changes". \emph{Proceedings of
    the Seventeenth Conference on Uncertainty in Artificial Intelligence},
    512--521.
}
\author{Robert Osazuwa Ness}
\seealso{\link{cpdag}, \link{score}.}
\keyword{graphs}
