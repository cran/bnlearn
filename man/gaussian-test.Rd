\name{gaussian.test}
\docType{data}
\alias{gaussian.test}
\title{Synthetic (continuous) data set to test learning algorithms}
\description{

  This a synthetic data set used as a test case in the \pkg{bnlearn} package.

}
\usage{
data(gaussian.test)
}
\format{

  The \code{gaussian.test} data set contains seven normal (Gaussian) variables.

}
\note{

  The R script to generate data from this network is available from
  \url{https://www.bnlearn.com/documentation/networks/}.

}
\examples{
# load the data.
data(gaussian.test)
# create and plot the network structure.
dag = model2network("[A][B][E][G][C|A:B][D|B][F|A:D:E:G]")
\dontrun{graphviz.plot(dag)}
}
\keyword{datasets}
