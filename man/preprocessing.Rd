\name{preprocess}
\alias{discretize}
\alias{dedup}
\title{Pre-process data to better learn Bayesian networks}
\description{

  Screen and transform the data to make them more suitable for structure and
  parameter learning.

}
\usage{
  # discretize continuous data into factors.
  discretize(data, method, breaks = 3, ordered = FALSE, ..., debug = FALSE)
  # screen continuous data for highly correlated pairs of variables.
  dedup(data, threshold, debug = FALSE)
}
\arguments{
  \item{data}{a data frame containing numeric columns (for \code{dedup}) or a
    combination of numeric or factor columns (for \code{}).}
  \item{threshold}{a numeric value between zero and one, the absolute
    correlation used a threshold in screening highly correlated pairs.}
  \item{method}{a character string, either \code{interval} for \emph{interval
    discretization}, \code{quantile} for \emph{quantile discretization}
    (the default) or \code{hartemink} for \emph{Hartemink's pairwise mutual
    information} method.}
  \item{breaks}{if \code{method} is set to \code{hartemink}, an integer number,
    the number of levels the variables are to be discretized into. Otherwise,
    a vector of integer numbers, one for each column of the data set, specifying
    the number of levels for each variable.}
  \item{ordered}{a boolean value. If \code{TRUE} the discretized variables are
    returned as ordered factors instead of unordered ones.}
  \item{\dots}{additional tuning parameters, see below.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
}
\details{

  \code{discretize} takes a data frame of continuous variables as its first
  argument and returns a secdond data frame of discrete variables, transformed
  using of three methods: \code{interval}, \code{quantile} or \code{hartemink}.

  \code{dedup} screens the data for pairs of highly correlated variables, and
   discards one in each pair.

}
\note{

  Hartemink's algorithm has been designed to deal with sets of homogeneous,
  continuous variables; this is the reason why they are initially transformed
  into discrete variables, all with the same number of levels (given by the
  \code{ibreaks} argument). Which of the other algorithms is used is specified
  by the \code{idisc} argument (\code{quantile} is the default). The
  implementation in \code{bnlearn} also handles sets of discrete variables
  with the same number of levels, which are treated as adjacent interval
  identifiers. This allows the user to perform the initial discretization
  with the algorithm of his choice, as long as all variables have the same
  number of levels in the end.

}
\value{

  \code{discretize} returns a data frame with the same structure (number of
  columns, column names, etc.) as \code{data}, containing the discretized
  variables.

  \code{dedup} returns a data frame with a subset of the columns of \code{data}.

}
\examples{
data(gaussian.test)
d = discretize(gaussian.test, method = 'hartemink', breaks = 4, ibreaks = 20)
plot(hc(d))
d2 = dedup(gaussian.test)
}
\references{

  Hartemink A (2001). \emph{Principled Computational Methods for the Validation
      and Discovery of Genetic Regulatory Networks}. Ph.D. thesis, School of
      Electrical Engineering and Computer Science, Massachusetts Institute of
      Technology.

}
\author{ Marco Scutari }
\keyword{manip}
