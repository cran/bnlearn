\name{rbn}
\alias{rbn}
\title{ Generate random data from a given Bayesian network }
\description{

  Generate random data from a given Bayesian network.

}
\usage{
rbn(x, n, data, debug = FALSE)
}
\arguments{
  \item{x}{an object of class \code{bn}.}
  \item{n}{non-negative integer giving the number of observations
       to generate.}
  \item{data}{a data frame, containing the data the Bayesian network
       was learned from.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output
       is printed; otherwise the function is completely silent.}
}
\value{

  A data frame with the same structure (column names and data types)
  of the \code{data} parameter.

}
\references{

  Korb K, Nicholson AE (2003). \emph{Bayesian Artificial Intelligence}.
      Chapman & Hall/CRC.

}
\examples{
\dontrun{
data(learning.test)
res = gs(learning.test)
res = set.arc(res, "A", "B")
par(mfrow = c(1,2))
plot(res)
sim = rbn(res, 500, learning.test)
plot(gs(sim))
}
}
\author{ Marco Scutari }
\keyword{utilities}
\keyword{multivariate}
