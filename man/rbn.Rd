\name{rbn}
\alias{rbn}
\alias{rbn.bn}
\alias{rbn.bn.fit}
\title{ Generate random data from a given Bayesian network }
\description{

  Generate random data from a given Bayesian network.

}
\usage{
\method{rbn}{bn}(x, n = 1, data, fit = "mle", \dots, debug = FALSE)
\method{rbn}{bn.fit}(x, n = 1, \dots, debug = FALSE)
}
\arguments{
  \item{x}{an object of class \code{bn} or \code{bn.fit}.}
  \item{n}{a positive integer giving the number of observations
       to generate.}
  \item{data}{a data frame containing the data the Bayesian network
       was learned from.}
  \item{fit}{a character string, the label of the method used to fit the
       parameters of the newtork. See \code{\link{bn.fit}} for details.}
  \item{...}{additional arguments for the parameter estimation prcoedure,
    see again \code{\link{bn.fit}} for details..}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output
       is printed; otherwise the function is completely silent.}
}
\value{

  A data frame with the same structure (column names and data types)
  of the \code{data} parameter.

}
\references{

  Korb K, Nicholson AE (2010). \emph{Bayesian Artificial Intelligence}.
      Chapman & Hall/CRC, 2nd edition.

}
\examples{
\dontrun{
data(learning.test)
res = gs(learning.test)
res = set.arc(res, "A", "B")
par(mfrow = c(1,2))
plot(res)
sim = rbn(res, 500, learning.test)
plot(gs(sim))}}
\author{ Marco Scutari }
\seealso{\code{\link{bn.boot}}, \code{\link{bn.cv}}.}
\keyword{utilities}
\keyword{multivariate}
