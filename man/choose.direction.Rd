\name{choose.direction}
\alias{choose.direction}
\title{ Try to infer the direction of an undirected arc }
\description{

  Check both possible directed arcs for existence,
  and choose the one with the lowest p-value.

}
\usage{

  choose.direction(x, arc, data, debug = FALSE)

}
\arguments{
  \item{x}{an object of class "bn".}
  \item{arc}{ a character string array of length 2, the labels of
       two nodes of the graph.}
  \item{data}{a data frame, containing the data the Bayesian network
       was learned from.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output 
       is printed; otherwise the function is completely silent.}
}
\value{

  \code{choose.direction} returns invisibly an updated copy of \code{x}.

}
\examples{

data(learning.test)

res = gs(learning.test)
# the arc between E and F has no direction
plot(res)
res = choose.direction(res, c("E", "F"), learning.test, debug = TRUE)
# * testing E - F for direction.
#   > testing E -> F with conditioning set '  '.
#     > p-value is 2.174867e-197 .
#   > testing F -> E with conditioning set ' B '.
#     > p-value is 3.935648e-245 .
#   @ removing E -> F .
plot(res, highlight = c("E", "F"))

}
\author{ Marco Scutari }
\keyword{htest}
