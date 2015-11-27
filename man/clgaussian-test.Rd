\name{clgaussian.test}
\docType{data}
\alias{clgaussian.test}
\title{Synthetic (mixed) data set to test learning algorithms}
\description{

  This a synthetic data set used as a test case in the \pkg{bnlearn} package.

}
\usage{
data(gaussian.test)
}
\format{

  The \code{clgaussian.test} data set contains one normal (Gaussian) variable,
  4 discrete variables and 3 conditional Gaussian variables.

}
\note{

  The R script to generate data from this network is shipped in the
  \file{network.scripts} directory of this package.

}
\examples{
# load the data and build the correct network from the model string.
data(clgaussian.test)
res = empty.graph(names(clgaussian.test))
modelstring(res) = "[A][B][C][H][D|A:H][F|B:C][E|B:D][G|A:D:E:F]"
plot(res)
}
\keyword{datasets}
