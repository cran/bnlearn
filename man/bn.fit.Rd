\name{bn.fit}
\alias{bn.fit}
\title{ Fit the parameters of a Bayesian network }
\description{

  Fit the parameters of a Bayesian network conditional
  on its structure.

}
\usage{
  bn.fit(x, data, debug = FALSE)
}
\arguments{
  \item{x}{an object of class \code{bn}.}
  \item{data}{a data frame, containing the variables in the model.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output
       is printed; otherwise the function is completely silent.}
}
\value{

  An object of class \code{bn.fit}. See \code{\link{bn.fit class}} for
  details.

}
\author{ Marco Scutari }
\keyword{multivariate}
