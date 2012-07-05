\name{bn.boot}
\alias{bn.boot}
\title{ Parametric and nonparametric bootstrap of Bayesian networks }
\description{

  Apply a user-specified function to the Bayesian network structures
  learned from bootstrap samples of the original data.

}
\usage{
bn.boot(data, statistic, R = 200, m = nrow(data),
  sim = "ordinary", algorithm, algorithm.args = list(),
  statistic.args = list(), cluster = NULL, debug = FALSE)
}
\arguments{
  \item{data}{a data frame containing the variables in the model.}
  \item{statistic}{a function or a character string (the name of a
    function) to be applied to each bootstrap replicate.}
  \item{R}{a positive integer, the number of bootstrap replicates.}
  \item{m}{a positive integer, the size of each bootstrap replicate.}
  \item{sim}{a character string indicating the type of simulation
    required. Possible values are \code{"ordinary"} (the default)
    and \code{"parametric"}.}
  \item{algorithm}{a character string, the learning algorithm to be
    applied to the bootstrap replicates. Possible values are \code{gs},
    \code{iamb}, \code{fast.iamb}, \code{inter.iamb}, \code{mmpc},
    \code{hc}, \code{tabu}, \code{mmhc} and \code{rsmax2}. See
    \code{\link{bnlearn-package}} and documentation of each algorithm
    for details.}
  \item{algorithm.args}{a list of extra arguments to be passed to
    the learning algorithm.}
  \item{statistic.args}{a list of extra arguments to be passed to
    the function specified by \code{statistic}.}
  \item{cluster}{an optional cluster object from package \pkg{snow}.
      See \code{\link{snow integration}} for details and a simple
      example.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output
       is printed; otherwise the function is completely silent.}
}
\details{

  The first argument of \code{statistic} is the \code{bn} object encoding
  the network structure learned from the bootstrap sample; the arguments
  specified in \code{statistics.args} are extracted from the list and
  passed to \code{statitstics} as the 2nd, 3rd, etc. arguments.

}
\value{

  A list containing the results of the calls to \code{statistic}.

}
\references{

  Friedman N, Goldszmidt M, Wyner A (1999). "Data Analysis with Bayesian
      Networks: A Bootstrap Approach". In "UAI '99: Proceedings of the
      15th Annual Conference on Uncertainty in Artificial Intelligence",
      pp. 196-20. Morgan Kaufmann.

}
\examples{
\dontrun{
data(learning.test)
bn.boot(data = learning.test, R = 2, m = 500, algorithm = "gs",
  statistic = arcs)
# [[1]]
#      from to
# <arcs for the first replicate>
#
# [[2]]
#      from to
# <arcs for the second replicate>
}
}
\author{ Marco Scutari }
\seealso{\code{\link{bn.cv}}, \code{\link{rbn}}. }
\keyword{nonparametric}
\keyword{multivariate}
