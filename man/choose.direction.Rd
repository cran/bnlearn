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
  \item{data}{a data frame, containing the variables in the model.}
  \item{debug}{a boolean value. If TRUE a lot of debugging output is
       printed; otherwise the function is completely silent.}
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
# > testing E -> F with conditioning set '  '.
#   > p-value is 2.041687e-44 .
# > testing F -> E with conditioning set '  '.
#   > p-value is 4.598853e-08 .
# @ removing F -> E .
plot(res, highlight = c("E", "F"))

}
\author{ Marco Scutari }
\keyword{htest}
