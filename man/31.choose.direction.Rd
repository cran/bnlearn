\name{choose.direction}
\alias{choose.direction}
\title{ Try to infer the direction of an undirected arc }
\description{

  Check both possible directed arcs for existence,
  and choose the one with the lowest p-value.

}
\usage{
  choose.direction(x, arc, data, criterion = NULL, ..., debug = FALSE)
}
\arguments{
  \item{x}{an object of class "bn".}
  \item{arc}{ a character string vector of length 2, the labels of
       two nodes of the graph.}
  \item{data}{a data frame, containing the data the Bayesian network
       was learned from.}
  \item{criterion}{the label of either a score or an independence test.
       See \code{\link{bnlearn-package}} for details.}
  \item{...}{additional tuning parameters for the network score.
      See \code{\link{score}} for details.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output 
       is printed; otherwise the function is completely silent.}
}
\value{

  \code{choose.direction} returns invisibly an updated copy of \code{x}.

}
\examples{
data(learning.test)
res = gs(learning.test)

## the arc A - B has no direction.
choose.direction(res, data = learning.test, arc = c("A", "B"), debug = TRUE)
# * testing A - B for direction.
#   > testing A -> B with conditioning set '  '.
#     > p-value is 0 .
#   > testing B -> A with conditioning set '  '.
#     > p-value is 0 .
#   @ nothing to do, same p-value.

## let's see score equivalence in action.
choose.direction(res, data = learning.test, criterion = "aic", 
  arc = c("A", "B"), debug = TRUE)
# * testing A - B for direction.
#     > delta between scores for nodes A B is 1166.914 .
#     > delta between scores for nodes B A is 1166.914 .
#   > initial score for node A is -5495.051 .
#   > initial score for node B is -4834.284 .
#   > score delta for arc A -> B is 1166.914 .
#   > score delta for arc B -> A is 1166.914 .
#   @ nothing to do, same score delta.

## arcs which introduce cycles are handled correctly.
choose.direction(res, data = learning.test, arc = c("A", "E"), debug = TRUE)
# * testing A - E for direction.
#   > testing A -> E with conditioning set ' B F '.
#     > p-value is 0.9818423 .
#   > testing E -> A with conditioning set '  '.
#    > p-value is 1.426725e-99 .
#   > adding E -> A creates cycles!.
#   > arc A -> E isn't good, either.
}
\author{ Marco Scutari }
\keyword{htest}
