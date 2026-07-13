\name{data preprocessing}
\alias{discretize}
\alias{dedup}
\title{Pre-process data to better learn Bayesian networks}
\description{

  Screen and transform the data to make them more suitable for structure and
  parameter learning.

}
\usage{
  # discretise continuous data into factors.
  discretize(data, method, breaks = 3, ordered = FALSE, ..., debug = FALSE)
  # screen continuous data for highly correlated pairs of variables.
  dedup(data, method, threshold, debug = FALSE)
}
\arguments{
  \item{data}{a data frame containing numeric columns (for \code{dedup()}) or a
    combination of numeric or factor columns (for \code{discretize()}).}
  \item{threshold}{a numeric value between zero and one, the absolute
    correlation used as a threshold in screening highly correlated pairs.}
  \item{method}{a character string, the label of the method used to preprocess
    the data. For \code{discretize()}, the possible values are \code{"interval"}
    for \emph{interval discretisation}, \code{"quantile"} for \emph{quantile
    discretisation} (the default) or \code{"hartemink"} for \emph{Hartemink's
    pairwise mutual information} method. For \code{dedup()}, the only possible
    value is \code{"cor"} for screening based on \emph{linear correlation}.}
  \item{breaks}{an integer number, the number of levels the variables will be
    discretised into; or a vector of integer numbers, one for each column of the
    data set, specifying the number of levels for each variable.}
  \item{ordered}{a boolean value. If \code{TRUE}, the discretised variables are
    returned as ordered factors instead of unordered ones.}
  \item{\dots}{additional tuning parameters, see below.}
  \item{debug}{a boolean value. If \code{TRUE}, a lot of debugging output is
    printed. Otherwise, the function is completely silent.}
}
\details{

  \code{discretize()} takes a data frame as its first argument and returns a
  second data frame of discrete variables, transformed using one of three
  methods: \code{"interval"}, \code{"quantile"} or \code{"hartemink"}. Discrete
  variables are left unchanged.

  The \code{"hartemink"} method has two additional tuning parameters:
  \itemize{

    \item \code{idisc}: the method used for the initial marginal discretisation
      of the variables, either \code{"interval"} or \code{"quantile"}.
    \item \code{ibreaks}: the number of levels the variables are initially
      discretised into, in the same format as in the \code{breaks} argument.

  }

  It is sometimes the case that the \code{"quantile"} method cannot discretise one
  or more variables in the data without generating zero-length intervals because
  the quantiles are not unique. If \code{method = "quantile"},
  \code{discretize()} will produce an error. If \code{method = "quantile"} and
  \code{idisc = "quantile"}, \code{discretize()} will try to lower the number of
  breaks set by the \code{ibreaks} argument until quantiles are distinct. If
  this is not possible without making \code{ibreaks} smaller than \code{breaks},
  \code{discretize()} will return an error.

  \code{dedup()} screens the data for pairs of highly correlated variables, and
  discards one in each pair.

  Both \code{discretize()} and \code{dedup()} accept data with missing values.

}
\value{

  \code{discretize()} returns a data frame with the same structure (number of
  columns, column names, etc.) as \code{data}, containing the discretised
  variables. The data frame has an attribute \code{"cutpoints"} containing the
  cutpoints used internally by \code{discretize()}, which are made available so
  that discretising another data set yields the same set of factors.

  \code{dedup()} returns a data frame with a subset of the columns of \code{data}.

}
\examples{
data(gaussian.test)
d = discretize(gaussian.test, method = 'hartemink', breaks = 4, ibreaks = 10)
plot(hc(d))
d2 = dedup(gaussian.test)
}
\references{

  Hartemink A (2001). \emph{Principled Computational Methods for the Validation
    and Discovery of Genetic Regulatory Networks}. Ph.D. thesis, School of
    Electrical Engineering and Computer Science, Massachusetts Institute of
    Technology.

}
\author{Marco Scutari}
\keyword{data preprocessing}
