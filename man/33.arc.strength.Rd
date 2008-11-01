\name{arc.strength}
\alias{arc.strength}
\title{ Measure the strength of the arcs }
\description{

  Strength of the probabilistic relations expressed by one or more
  arcs of the Bayesian network.

}
\usage{
  arc.strength(x, data, criterion = NULL, ..., debug = FALSE)
}
\arguments{
  \item{x}{an object of class \code{bn}.}
  \item{data}{a data frame, containing the data the Bayesian network
       was learned from.}
  \item{criterion}{the label of either a score or an independence test.
       See \code{\link{bnlearn-package}} for details.}
  \item{\dots}{additional tuning parameters for the network score.
      See \code{\link{score}} for details.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output
       is printed; otherwise the function is completely silent.}
}
\value{

  \code{arc.strength} returns a data frame with three columns: \code{from},
  \code{to} and \code{strength}.
  If \code{criterion} is a conditional indipendence test the strength
  is a p-value (so the lower the value, the stronger the relationship);
  otherwise it's the score gain/loss which would be caused by the
  arc's removal.

}
\examples{
data(learning.test)
res = gs(learning.test)
res = set.arc(res, "A", "B")
arc.strength(res, learning.test)
#   from to      strength
# 1    A  B  0.000000e+00
# 2    A  D  0.000000e+00
# 3    B  E 1.024198e-320
# 4    C  D  0.000000e+00
# 5    F  E 3.935648e-245
arc.strength(res, learning.test, criterion = "aic")
#   from to   strength
# 1    A  B -1166.9139
# 2    A  D -1978.0531
# 3    B  E  -746.8954
# 4    C  D  -862.8637
# 5    F  E  -568.7816
res = set.arc(res, "B", "A")
# A -> B and B -> A have the same strength because they
# are score equivalent.
arc.strength(res, learning.test, criterion = "aic")
#   from to   strength
# 1    A  D -1978.0531
# 2    B  E  -746.8954
# 3    C  D  -862.8637
# 4    F  E  -568.7816
# 5    B  A -1166.9139
}
\seealso{\code{\link{choose.direction}}, \code{\link{score}}. }
\author{ Marco Scutari }
\keyword{htest}
