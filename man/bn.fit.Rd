\name{bn.fit}
\alias{bn.fit}
\title{ Fit the parameters of a Bayesian network }
\description{

  Fit the parameters of a Bayesian network conditional
  on its structure.

}
\usage{
bn.fit(x, data, method = "mle", \dots, debug = FALSE)
}
\arguments{
  \item{x}{an object of class \code{bn}.}
  \item{data}{a data frame containing the variables in the model.}
  \item{method}{a character string, either \code{mle} for \emph{Maximum
    Likelihood parameter estimation} or \code{bayes} for \emph{Bayesian
    parameter estimation} (currently implemented only for discrete data).}
  \item{...}{additional arguments for the parameter estimation prcoedure,
    see below.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output
    is printed; otherwise the function is completely silent.}
}
\note{

  Due to the way Bayesian networks are defined it's possible to
  estimate their parameters only if the network structure is
  completely directed (i.e. there are no undirected arcs). See
  \code{\link{set.arc}} and \code{\link{pdag2dag}} for two ways
  of manually setting the direction of one or more arcs.

  The only supported additional parameter is the imaginary sample size
  (\code{iss}) for the Dirichlet posterior distribution of discrete
  networks (see \code{\link{score}} for details).

}
\value{

  An object of class \code{bn.fit}. See \code{\link{bn.fit class}} for
  details.

}
\seealso{
  \code{\link{bn.fit utilities}}, \code{\link{bn.fit plots}}.
}
\examples{
data(learning.test)

# learn the network structure.
res = gs(learning.test)
# set the direction of the only undirected arc, A - B.
res = set.arc(res, "A", "B")
# estimate the parameters of the Bayesian network.
fitted = bn.fit(res, learning.test)
}
\author{ Marco Scutari }
\keyword{multivariate}