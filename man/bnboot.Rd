\name{bnboot}
\alias{bnboot}
\title{ Paramentric and nonparametric bootstrap of Bayesian networks }
\description{

  Apply a user-specified function to Bayesian networks learned
  from bootstrap relicates of the original data.

}
\usage{
  bnboot(data, statistic, R = 200, m = round(2/3 * nrow(data)),
    sim = "ordinary", algorithm, algorithm.args = list(),
    statistic.args = list(), debug = FALSE)
}
\arguments{
  \item{data}{a data frame, containing the variables in the model.}
  \item{statistic}{a function or a character string (the name of a
    function) to be applied to each bootstrap replicate.}
  \item{R}{a positive integer, the number of bootstrap replicates.}
  \item{m}{a positive integer, the size of each bootstrap replicate.}
  \item{sim}{a character string indicating the type of simulation
    required. Possible values are \code{"ordinary"} (the default)
    and \code{"parametric"}.}
  \item{algorithm}{a character string, the learning algorithm to be
    applied to the bootstrap replicates. Possible values are \code{gs},
    \code{iamb}, \code{fast.iamb}, \code{inter.iamb}, \code{mmpc}
    and \code{hc}. See \code{\link{bnlearn-package}} and
    documentation of each algorithm for details.}
  \item{algorithm.args}{a list of extra arguments to be passed to
    the learning algorithm.}
  \item{statistic.args}{a list of extra arguments to be passed to
    the function specified by \code{statistic}.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output
       is printed; otherwise the function is completely silent.}
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
bnboot(data = learning.test, R = 2, m = 500, algorithm = "gs",
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
\seealso{\code{\link{constraint-based algorithms}}, \code{\link{mmpc}}, 
  \code{\link{hc}}. }
\keyword{nonparametric}
\keyword{multivariate}
