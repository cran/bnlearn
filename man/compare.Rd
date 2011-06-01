\name{compare}
\alias{compare}
\alias{all.equal.bn}
\alias{shd}
\alias{hamming}
\title{ Compare two different Bayesian networks }
\description{

 Compare two different Bayesian networks; compute the Structural
 Hamming Distance (SHD) between them or the Hamming distance
 between their skeletons.

}
\usage{
compare(target, current, debug = FALSE)
\method{all.equal}{bn}(target, current, ...)

shd(learned, true, debug = FALSE)
hamming(learned, true, debug = FALSE)
}
\arguments{
  \item{target, learned}{an object of class \code{bn}.}
  \item{current, true}{another object of class \code{bn}.}
  \item{\dots}{extra arguments from the generic method (currently ignored).}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output
       is printed; otherwise the function is completely silent.}
}
\note{

  The main difference between \code{compare} and \code{all.equal} is
  that the former is quite slow but provides a detailed description
  of the differences between \code{target} and \code{current} (when
  \code{debug} is set to \code{TRUE}), while the latter is much faster
  but provides only generic diagnostic messages.

}
\value{

  \code{compare} returns a boolean value (\code{TRUE} if the objects
  describe the same network structure, \code{FALSE} otherwise).
  \code{all.equal} returns either \code{TRUE} or a character string
  describing the differences between \code{target} and \code{current}.
  \code{shd} returns a non-negative integer number.

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
\keyword{graphs}
