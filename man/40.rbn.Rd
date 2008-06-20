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
\note{

  Only discrete Bayesian networks are supported; the simulation
  uses an implementation of the Logic Sampling (LS) algorithm
  (Korb and Nicholson, 2004).

  The execution time scales linearly with the number of generated data.

  Node ordering uses a recursive approach; if there are lots of nodes
  in the network, you will need to increase R's recursion limit.  See
  the \code{expression} argument to the \code{options} command for
  details on how to do this.

}
\value{

  A data frame with the same structure (column names and data types)
  of the \code{data} parameter.

}
\references{

  K. Korb and A. Nicholson. Bayesian artificial intelligence. Chapman and
   Hall, 2004.

}
\examples{
\dontrun{
library(bnlearn)
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
