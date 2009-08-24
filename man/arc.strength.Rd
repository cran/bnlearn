\name{arc.strength}
\alias{arc.strength}
\title{ Measure the strength of the arcs }
\description{

  Strength of the probabilistic relations expressed by the
  arcs of the Bayesian network.

}
\usage{
  arc.strength(x, data, criterion = NULL, ..., debug = FALSE)
}
\arguments{
  \item{x}{an object of class \code{bn}.}
  \item{data}{a data frame, containing the data the Bayesian network
       was learned from.}
  \item{criterion}{a character string, the label of a score
       function, the label of an independence test or \code{bootstrap}.
       See \code{\link{bnlearn-package}} for details on the first two
       possibilities.}
  \item{\dots}{additional tuning parameters for the network score (if
      \code{criterion} is the label of a score function, see
      \code{\link{score}} for details), the conditional independence test
      (currently the only one is \code{B}, the number of permutations)
      or the bootstrap simulation (if \code{criterion} is set to
      \code{bootstrap}, see \code{\link{boot.strength}} for details).}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output
       is printed; otherwise the function is completely silent.}
}
\details{

  If \code{criterion} is a conditional independence test the strength
  is a p-value (so the lower the value, the stronger the relationship).
  The only possible additional parameter is \code{B}, the number of
  permutations to be generated for a permutation test.

  If \code{criterion} is a score function it's the score gain/loss which
  would be caused by the arc's removal. There may be additional parameters
  depending on the choice of the score, see \code{\link{score}} for details.

  If \code{criterion} is \code{bootstrap} the strength is computed as
  in \code{\link{boot.strength}}. The additional parameters are \code{R},
  \code{m}, \code{algorithm} and \code{algorithm.args}; if the latter two
  are not specified, the values stored in \code{x} are used.

}
\value{

  \code{arc.strength} returns an object of class \code{bn.strength}.
  See \code{\link{bn.strength class}} for details.

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
\dontrun{
arc.strength(res, data = learning.test, criterion = "bootstrap",
  R = 200, algorithm.args = list(alpha = 0.10))
#   from to strength
# 1    A  B        1
# 2    A  D        1
# 3    B  E        1
# 4    C  D        1
# 5    F  E        1
}
}
\seealso{\code{\link{choose.direction}}, \code{\link{score}}. }
\author{ Marco Scutari }
\keyword{htest}
