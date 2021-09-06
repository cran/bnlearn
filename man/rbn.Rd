\name{rbn}
\alias{rbn}
\alias{rbn.bn}
\alias{rbn.bn.fit}
\title{Simulate random samples from a given Bayesian network}
\description{

  Simulate random samples from a given Bayesian network.

}
\usage{
\method{rbn}{bn}(x, n = 1, data, fit = "mle", \dots, debug = FALSE)
\method{rbn}{bn.fit}(x, n = 1, \dots, debug = FALSE)
}
\arguments{
  \item{x}{an object of class \code{bn} or \code{bn.fit}.}
  \item{n}{a positive integer giving the number of observations to generate.}
  \item{data}{a data frame containing the data the Bayesian network was learned
    from.}
  \item{fit}{a character string, the label of the method used to fit the
    parameters of the newtork. See \code{\link{bn.fit}} for details.}
  \item{...}{additional arguments for the parameter estimation prcoedure, see
    again \code{\link{bn.fit}} for details.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
}
\details{

  \code{rbn()} implements forward/logic sampling: values for the root nodes are
  sampled from their (unconditional) distribution, then those of their children
  conditional on the respective parent sets. This is done iteratively until
  values have been sampled for all nodes.

  If \code{x} contains \code{NA} parameter estimates (because of unobserved
  discrete parents configurations in the data the parameters were learned from),
  \code{rbn} will produce samples that contain \code{NA}s when those parents
  configurations appear in the simulated samples. See \code{\link{bn.fit}} for
  details on how to make sure \code{bn.fit} objects contain no \code{NA}
  parameter estimates.

}
\value{

  A data frame with the same structure (column names and data types) of the
  \code{data} argument (if \code{x} is an object of class \code{bn}) or with
  the same structure as the data originally used to to fit the parameters of
  the Bayesian network (if \code{x} is an object of class \code{bn.fit}).

}
\references{

  Korb K, Nicholson AE (2010). \emph{Bayesian Artificial Intelligence}.
    Chapman & Hall/CRC, 2nd edition.

}
\examples{
\dontrun{
data(learning.test)
res = hc(learning.test)
par(mfrow = c(1, 2))
plot(res)
sim = rbn(res, 500, learning.test)
plot(hc(sim))}}
\author{Marco Scutari}
\seealso{\code{\link{bn.boot}}, \code{\link{bn.cv}}.}
\keyword{inference}
\keyword{simulation}
