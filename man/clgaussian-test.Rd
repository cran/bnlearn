\name{clgaussian.test}
\docType{data}
\alias{clgaussian.test}
\title{Synthetic (mixed) data set to test learning algorithms}
\description{

  This a synthetic data set used as a test case in the \pkg{bnlearn} package.

}
\usage{
data(clgaussian.test)
}
\format{

  The \code{clgaussian.test} data set contains one normal (Gaussian) variable,
  4 discrete variables and 3 conditional Gaussian variables.

}
\note{

  The R script to generate data from this network is available from
  \url{https://www.bnlearn.com/documentation/networks/}.

}
\examples{
# load the data.
data(clgaussian.test)
# create and plot the network structure.
dag = model2network("[A][B][C][H][D|A:H][F|B:C][E|B:D][G|A:D:E:F]")
\dontrun{graphviz.plot(dag)}
}
\keyword{datasets}
