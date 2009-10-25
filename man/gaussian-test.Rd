\name{gaussian.test}
\docType{data}
\alias{gaussian.test}
\title{ Synthetic (continuous) data set to test learning algorithms }
\description{

  This a synthetic data set used as a test case in the
  \pkg{bnlearn} package.

}
\usage{
data(gaussian.test)
}
\format{

  The \code{gaussian.test} data set contains the seven normal
  (Gaussian) variables.

}
\note{

  The R script to generate data from this network is shipped in the
  \file{network.scripts} directory of this package.

}
\examples{
# load the data and build the correct network from the model string.
data(gaussian.test)
res = empty.graph(names(gaussian.test))
modelstring(res) = "[A][B][E][G][C|A:B][D|B][F|A:D:E:G]"
plot(res)
}
\keyword{datasets}
