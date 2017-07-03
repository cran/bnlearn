\name{structural.em}
\alias{em-based algorithms}
\alias{structural.em}
\title{Structure learning from missing data}
\description{

  Learn the structure of a Bayesian network from a data set containing missing
  values using Structural EM.

}
\usage{
structural.em(x, maximize = "hc", maximize.args = list(), fit = "mle", 
    fit.args = list(), impute, impute.args = list(), return.all = FALSE,
    max.iter = 5, debug = FALSE)
}
\arguments{
  \item{x}{a data frame containing the variables in the model.}
  \item{maximize}{a character string, the score-based algorithm to be used in
    the \dQuote{maximization} step. See \code{\link{bnlearn-package}} for
    details.}
  \item{maximize.args}{a list of arguments to be passed to the algorithm
    specified by \code{maximize}, such as \code{restart} for hill-climbing or
    \code{tabu} for tabu search.}
  \item{fit}{a character string, the parameter learning method to be used in
    the \dQuote{maximization} step. See \code{\link{bn.fit}} for details.}
  \item{fit.args}{a list of arguments to be passed to the parameter learning
    method specified by \code{fit}.}
  \item{impute}{a character string, the imputation method to be used in the
    \dQuote{expectation} step. See \code{\link{impute}} for details.}
  \item{impute.args}{a list of arguments to be passed to the imputation method
    specified by \code{impute}.}
  \item{return.all}{a boolean value. See below for details.}
  \item{max.iter}{an integer, the maximum number of iterations.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
}
\value{

  If \code{return.all} is \code{FALSE}, \code{structural.em} returns an object
  of class \code{bn}. (See \code{\link{bn-class}} for details.)

  If \code{return.all} is \code{TRUE}, \code{structural.em} returns a list with
  three elements named \code{dag} (an object of class \code{bn}), \code{imputed}
  (a data frame containing the imputed data from the last iteration) and
  \code{fitted} (an object of class \code{bn.fit}, again from the last
  iteration; see \code{\link{bn.fit-class}} for details).

}
\references{

  Friedman N (1997). "Learning Belief Networks in the Presence of Missing Values
    and Hidden Variables". In "Proceedings of the 14th International Conference
    on Machine Learning (ICML)", pp. 125-133.

}
\author{Marco Scutari}
\seealso{\link{score-based algorithms}, \link{bn.fit}, \link{impute}.}
\keyword{structure learning}
\keyword{missing data}
