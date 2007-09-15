\name{learning.test}
\docType{data}
\alias{learning.test}
\title{ Synthetic dataset to test learning algorithms }
\description{

  This a synthetic dataset used as a test case in the
  bnlearn package.

}
\usage{

  learning.test

}
\examples{

data(learning.test)
str(learning.test)
# 'data.frame':   5000 obs. of  6 variables:
#  $ A: Factor w/ 3 levels "a","b","c": 2 2 1 1 1 3 3 2 2 2 ...
#  $ B: Factor w/ 3 levels "a","b","c": 3 1 1 1 1 3 3 2 2 1 ...
#  $ C: Factor w/ 3 levels "a","b","c": 2 3 1 1 2 1 2 1 2 2 ...
#  $ D: Factor w/ 3 levels "a","b","c": 1 1 1 1 3 3 3 2 1 1 ...
#  $ E: Factor w/ 3 levels "a","b","c": 2 2 1 2 1 3 3 2 3 1 ...
#  $ F: Factor w/ 2 levels "a","b": 2 2 1 2 1 1 1 2 1 1 ...

gs(learning.test)

}
\keyword{datasets}
