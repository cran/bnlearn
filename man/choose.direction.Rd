\name{choose.direction}
\alias{choose.direction}
\title{ Try to infer the direction of an undirected arc }
\description{

  Check both possible directed arcs for existence,
  and choose the one with the lowest p-value, the
  highest score or the highest bootstrap probability.

}
\usage{
choose.direction(x, arc, data, criterion = NULL, ...,
  debug = FALSE)
}
\arguments{
  \item{x}{an object of class \code{bn}.}
  \item{arc}{ a character string vector of length 2, the labels of
       two nodes of the graph.}
  \item{data}{a data frame containing the data the Bayesian network
       was learned from.}
  \item{criterion}{a character string, the label of a score
       function, the label of an independence test or \code{bootstrap}.
       See \code{\link{bnlearn-package}} for details on the first two
       possibilities.}
  \item{\dots}{additional tuning parameters for the network score.
      See \code{\link{score}} for details.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output
       is printed; otherwise the function is completely silent.}
}
\details{

  If \code{criterion} is \code{bootstrap}, \code{choose.directions} accepts
  the same arguments as \code{boot.strength}: \code{R} (the number of 
  bootstrap replicates), \code{m} (the bootstrap sample size), \code{algorithm}
  (the structure learning algorithm), \code{algorithm.args} (the arguments to
  pass to the structure learning algorithm) and \code{cpdag} (whether to
  transform the network structure to the CPDAG representation of the equivalence
  class it belongs to).

}
\value{

  \code{choose.direction} returns invisibly an updated copy of \code{x}.

}
\examples{
data(learning.test)
res = gs(learning.test)

## the arc A - B has no direction.
choose.direction(res, learning.test, arc = c("A", "B"), debug = TRUE)
# * testing A - B for direction.
#   > testing A -> B with conditioning set '  '.
#     > p-value is 0 .
#   > testing B -> A with conditioning set '  '.
#     > p-value is 0 .
#   @ nothing to do, same p-value.

## let's see score equivalence in action.
choose.direction(res, learning.test, criterion = "aic",
  arc = c("A", "B"), debug = TRUE)
# * testing A - B for direction.
#   > initial score for node A is -5495.051 .
#   > initial score for node B is -4834.284 .
#   > score delta for arc A -> B is 1166.914 .
#   > score delta for arc B -> A is 1166.914 .
#   @ nothing to do, same score delta.

## arcs which introduce cycles are handled correctly.
res = set.arc(res, "A", "B")
# now A -> B -> E -> A is a cycle.
choose.direction(res, learning.test, arc = c("E", "A"), debug = TRUE)
# * testing E - A for direction.
#   > testing E -> A with conditioning set '  '.
#     > p-value is 1.426725e-99 .
#   > testing A -> E with conditioning set ' B F '.
#     > p-value is 0.9818423 .
#   > adding E -> A creates cycles!.
#   > arc A -> E isn't good, either.

\dontrun{
choose.direction(res, learning.test, arc = c("D", "E"), criterion = "bootstrap",
  R = 100, algorithm = "iamb", algorithm.args = list(test = "x2"), cpdag = TRUE,
  debug = TRUE)
# * testing D - E for direction.
#   > testing D -> E 
#     > bootstrap probability of an arc between D and E is 0.36 .
#     > direction confidence for arc D -> E is 0.3333333 .
#     > direction confidence for arc E -> D is 0.6666667 .
#   @ nothing to do, bootstrap probability is less than 0.50.
}}
\seealso{\code{\link{score}}, \code{\link{arc.strength}}. }
\author{ Marco Scutari }
\keyword{htest}
