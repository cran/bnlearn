\name{compare}
\alias{compare}
\title{ Compare two different Bayesian networks }
\description{

 Compare two different Bayesian networks.

}
\usage{

  compare(r1, r2, debug = FALSE)

}
\arguments{
  \item{r1}{an object of class "bn".}
  \item{r2}{another object of class "bn".}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output 
       is printed; otherwise the function is completely silent.}
}
\value{

  A boolean value: \code{TRUE} if the objects describe the same
  network structure, \code{FALSE} otherwise.

}
\examples{

data(learning.test)

res = gs(learning.test)
# the arc between A and B has no direction
plot(res, highlight = c("A", "B"))
res2 = set.arc(res, "A", "B")
compare(res, res2, debug = TRUE)
# * arcs in r1 not present in r2:
# [1] "B -> A"
# * arcs in r2 not present in r1:
# character(0)
# [1] FALSE

}
\author{ Marco Scutari }
\keyword{utilities}
