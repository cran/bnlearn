\name{multivariate normal distribution}
\alias{mvnorm2gbn}
\alias{gbn2mvnorm}
\title{Gaussian Bayesian networks and multivariate normals}
\description{

  Convert a Gaussian Bayesian network into the multivariate normal distribution
  that is its global distribution, and vice versa.

}
\usage{
gbn2mvnorm(fitted)
mvnorm2gbn(dag, mu, sigma)
}
\arguments{
  \item{fitted}{an object of class \code{bn.fit}.}
  \item{dag}{an object of class \code{bn}, the structure of the network that
    will be returned.}
  \item{mu}{a numeric vector, the expectation of the multivariate normal.}
  \item{sigma}{a square numeric matrix, the covariance matrix of the
    multivariate normal.}
}
\value{

  \code{gbn2mvnorm()} returns a list with elements \code{"mu"} (the vector of
  expectations) and \code{"sigma"} (the covariance matrix).

  \code{mvnorm2gbn()} returns an object of class \code{bn.fit}.

}
\references{

  Pourahmadi M (2011). "Covariance Estimation: The GLM and Regularization
    Perspectives". \emph{Statistical Science}, 26(3), 369--387.

}
\examples{
data(gaussian.test)
dag = model2network("[A][B][E][G][C|A:B][D|B][F|A:D:E:G]")
bn = bn.fit(dag, gaussian.test)
mvn = gbn2mvnorm(bn)
bn2 = mvnorm2gbn(dag, mu = mvn$mu, sigma = mvn$sigma)
all.equal(bn, bn2)
}
\author{Marco Scutari}
\seealso{\link{bn.fit}.}
\keyword{parameter learning}
