\name{learning.test}
\docType{data}
\alias{learning.test}
\title{ Synthetic (discrete) data set to test learning algorithms }
\description{

  This a synthetic data set used as a test case in the
  \pkg{bnlearn} package.

}
\usage{
data(learning.test)
}
\format{

  The \code{learning.test} data set contains the following variables:
  \itemize{

    \item \code{A}, a three-level factor with levels \code{a},
            \code{b} and \code{c}.
    \item \code{B}, a three-level factor with levels \code{a},
            \code{b} and \code{c}.
    \item \code{C}, a three-level factor with levels \code{a},
            \code{b} and \code{c}.
    \item \code{D}, a three-level factor with levels \code{a},
            \code{b} and \code{c}.
    \item \code{E}, a three-level factor with levels \code{a},
            \code{b} and \code{c}.
    \item \code{F}, a two-level factor with levels \code{a}
            and \code{b}.

  }

}
\note{

  The R script to generate data from this network is shipped in the
  \file{network.scripts} directory of this package.

}
\examples{
# load the data and build the correct network from the model string.
data(learning.test)
res = empty.graph(names(learning.test))
modelstring(res) = "[A][C][F][B|A][D|A:C][E|B:F]"
plot(res)
}
\keyword{datasets}
