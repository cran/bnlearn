\name{compare}
\alias{compare}
\alias{shd}
\title{ Compare two different Bayesian networks }
\description{

 Compare two different Bayesian networks or compute the Structural
 Hamming Distance (SHD) between them.

}
\usage{
compare(r1, r2, debug = FALSE)
shd(learned, true, debug = FALSE)
}
\arguments{
  \item{r1, learned}{an object of class \code{bn}.}
  \item{r2, true}{another object of class \code{bn}.}
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
## the arc between A and B has no direction
plot(res, highlight = c("A", "B"))
res2 = set.arc(res, "A", "B")
compare(res, res2, debug = TRUE)
# * children of A in r1 not present in r2:
# character(0)
# * children of A in r2 not present in r1:
# [1] "B"
# * parents of B in r1 not present in r2:
# character(0)
# * parents of B in r2 not present in r1:
# [1] "A"
# * directed arcs in r1 not present in r2:
# character(0)
# * directed arcs in r2 not present in r1:
# [1] "A -> B"
# * undirected arcs in r1 not present in r2:
# [1] "A - B" "B - A"
# * undirected arcs in r2 not present in r1:
# character(0)
# [1] FALSE
e1 = model2network("[A][B][C|A:B][D|B][E|C][F|A:E]")
e2 = model2network("[A][B][C|A:B][D|B][E|C:F][F|A]")
shd(e2, e1, debug = TRUE)
# * arcs between A and F do not match.
# * arcs between E and F do not match.
# [1] 2

}
\references{

  Tsamardinos I, Brown LE, Aliferis CF (2006). "The Max-Min Hill-Climbing
      Bayesian Network Structure Learning Algorithm". \emph{Machine Learning},
      \strong{65}(1), 31-78.

}
\author{ Marco Scutari }
\keyword{utilities}
