\name{alpha.star}
\alias{alpha.star}
\title{Estimate the Optimal Imaginary Sample Size for BDe(u)}
\description{

  Estimate the optimal value of the imaginary sample size for the BDe score,
  assuming a uniform prior and given a network structure and a data set.

}
\usage{
  alpha.star(x, data, debug = FALSE)
}
\arguments{
  \item{x}{an object of class \code{bn} (for \code{bn.fit} and \code{custom.fit})
    or an object of class \code{bn.fit} (for \code{bn.net}).}
  \item{data}{a data frame containing the variables in the model.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
}
\value{

  \code{alpha.star} returns a positive number, the optimal imaginary sample
  size value.

}
\examples{
data(learning.test)
dag = hc(learning.test, score = "bic") 

for (i in 1:3) {

  a = alpha.star(dag, learning.test)
  dag = hc(learning.test, score = "bde", iss = a)

}#FOR
}
\references{

  Steck H (2008). "Learning the Bayesian Network Structure: Dirichlet Prior
      versus Data". In "Proceedings of the 24th Conference
      on Uncertainty in Artificial Intelligence (UAI2008)", pp. 511-518.

}
\author{ Marco Scutari }
\keyword{multivariate}
\keyword{models}
