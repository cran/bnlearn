\name{configs}
\alias{configs}
\title{Construct configurations of discrete variables}
\description{

  Create configurations of discrete variables, which can be used in modelling
  conditional probability tables.

}
\usage{
  configs(data, all = TRUE)
}
\arguments{
  \item{data}{a data frame containing factor columns.}
  \item{all}{a boolean value. If \code{TRUE} all configuration are included as
   levels in the return value;  otherwise only configurations which are actually
   observed are considered.}
}
\value{

  A factor with one element for each row of \code{data}, and levels as
  specified by \code{all}.

}
\examples{
data(learning.test)
configs(learning.test, all = TRUE)
configs(learning.test, all = FALSE)
}
\author{Marco Scutari}
\keyword{data preprocessing}
\keyword{convenience functions}
