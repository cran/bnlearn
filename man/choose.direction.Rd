\name{choose.direction}
\alias{choose.direction}
\title{Try to infer the direction of an undirected arc}
\description{

  Check both possible directed arcs for existence, and choose the one with the
  lowest p-value, the highest score or the highest bootstrap probability.

}
\usage{
choose.direction(x, arc, data, criterion = NULL, ..., debug = FALSE)
}
\arguments{
  \item{x}{an object of class \code{bn}.}
  \item{arc}{ a character string vector of length 2, the labels of two nodes of
    the graph.}
  \item{data}{a data frame containing the data the Bayesian network was learned
    from.}
  \item{criterion}{a character string, the label of a score function, the label
    of an independence test or \code{bootstrap}. See \code{\link{bnlearn-package}}
    for details on the first two possibilities.}
  \item{\dots}{additional tuning parameters for the network score. See
    \code{\link{score}} for details.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
}
\details{

  If \code{criterion} is \code{bootstrap}, \code{choose.directions} accepts the
  same arguments as \code{boot.strength}: \code{R} (the number of bootstrap
  replicates), \code{m} (the bootstrap sample size), \code{algorithm} (the
  structure learning algorithm), \code{algorithm.args} (the arguments to pass
  to the structure learning algorithm) and \code{cpdag} (whether to transform
  the network structure to the CPDAG representation of the equivalence class it
  belongs to).

}
\value{

  \code{choose.direction} returns invisibly an updated copy of \code{x}.

}
\examples{
data(learning.test)
res = gs(learning.test)

## the arc A - B has no direction.
choose.direction(res, learning.test, arc = c("A", "B"), debug = TRUE)

## let's see score equivalence in action.
choose.direction(res, learning.test, criterion = "aic",
  arc = c("A", "B"), debug = TRUE)

## arcs which introduce cycles are handled correctly.
res = set.arc(res, "A", "B")
# now A -> B -> E -> A is a cycle.
choose.direction(res, learning.test, arc = c("E", "A"), debug = TRUE)

\dontrun{
choose.direction(res, learning.test, arc = c("D", "E"), criterion = "bootstrap",
  R = 100, algorithm = "iamb", algorithm.args = list(test = "x2"), cpdag = TRUE,
  debug = TRUE)
}}
\seealso{\code{\link{score}}, \code{\link{arc.strength}}.}
\author{ Marco Scutari }
\keyword{htest}
