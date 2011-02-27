\name{discretize}
\alias{discretize}
\title{ Discretize data to learn discrete Bayesian networks }
\description{

  Discretize data to learn discrete Bayesian networks.

}
\usage{
  discretize(x, method, breaks = 3, ..., debug = FALSE)
}
\arguments{
  \item{x}{a data frame containing with numeric and factor columns.}  
  \item{method}{a character string, either \code{interval} for \emph{interval
    discretization}, \code{quantile} for \emph{quantile discretization} or
    \code{hartemink} for \emph{Hartemink's pairwise mutual information} method.}
  \item{breaks}{if \code{method} is set to \code{hartemink}, an integer number,
    the number of levels the variables are to be discretized into. Otherwise,
    a vector of integer numbers, one for each column of the data set, specifying
    the number of levels for each variable.}
  \item{\dots}{additional tuning parameters, see below.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output
      is printed; otherwise the function is completely silent.}
}
\value{

  \code{discretize} returns a data frame with the same structure (number
    of columns, column names, etc.) as \code{x}, containing the discretized
    variables.

}
\examples{
data(gaussian.test)
d = discretize(gaussian.test, method = 'hartemink', breaks = 4, initial.breaks = 20)
plot(hc(d))
}
\references{

  Hartemink A (2001). \emph{Principled Computational Methods for the Validation
      and Discovery of Genetic Regulatory Networks}. Ph.D. thesis, School of
      Electrical Engineering and Computer Science, Massachusetts Institute of
      Technology.

}
\author{ Marco Scutari }
\keyword{manip}
