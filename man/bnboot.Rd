\name{bn.boot}
\alias{bn.boot}
\title{Nonparametric bootstrap of Bayesian networks}
\description{

  Apply a user-specified function to the Bayesian network structures learned
  from bootstrap samples of the original data.

}
\usage{
bn.boot(data, statistic, R = 200, m = nrow(data), algorithm,
  algorithm.args = list(), statistic.args = list(), cluster,
  debug = FALSE)
}
\arguments{
  \item{data}{a data frame containing the variables in the model.}
  \item{statistic}{a function or a character string (the name of a function)
    to be applied to each bootstrap replicate.}
  \item{R}{a positive integer, the number of bootstrap replicates.}
  \item{m}{a positive integer, the size of each bootstrap replicate.}
  \item{algorithm}{a character string, the learning algorithm to be applied
    to the bootstrap replicates. See \code{\link{structure learning}} and the
    documentation of each algorithm for details.}
  \item{algorithm.args}{a list of extra arguments to be passed to the learning
    algorithm.}
  \item{statistic.args}{a list of extra arguments to be passed to the function
    specified by \code{statistic}.}
  \item{cluster}{an optional cluster object from package \pkg{parallel}.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
}
\details{

  The first argument of \code{statistic} is the \code{bn} object encoding the
  network structure learned from the bootstrap sample; the arguments specified
  in \code{statistics.args} are extracted from the list and passed to
  \code{statistics} as the 2nd, 3rd, etc. arguments.

}
\value{

  A list containing the results of the calls to \code{statistic}.

}
\references{

  Friedman N, Goldszmidt M, Wyner A (1999). "Data Analysis with Bayesian
	Networks: A Bootstrap Approach". \emph{Proceedings of the 15th Annual
	Conference on Uncertainty in Artificial Intelligence}, 196--201.

}
\examples{
\dontrun{
data(learning.test)
bn.boot(data = learning.test, R = 2, m = 500, algorithm = "gs",
  statistic = arcs)
}
}
\author{Marco Scutari}
\seealso{\code{\link{bn.cv}}, \code{\link{rbn}}.}
\keyword{resampling}
