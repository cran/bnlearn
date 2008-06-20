\name{gaussian.test}
\docType{data}
\alias{gaussian.test}
\title{ Synthetic (continuous) dataset to test learning algorithms }
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
\examples{

## The modelstring() of this data set is:
# [A][B][E][G][C|A:B:E][D|B][F|A:D:E:G]

# these are the R commands used to generate this data set.
\dontrun{
a = data.frame(A = rnorm(5000, 1, 1),
               B = rnorm(5000, 2, 3),
               C = rep(0, 5000),
               D = rep(0, 5000),
               E = rnorm(5000, 3.5, 2),
               F = rep(0, 5000),
               G = rnorm(5000, 3.5, 2))
a$C = 2 * (a$A + a$B) + rnorm(5000, 2, 0.5)
a$D = 8 * (a$B) + rnorm(5000, 6, 0.33)
a$F = 2 * a$A + a$D + a$E + 1.5 * a$G + rnorm(5000, 0, 1)
}

}
\keyword{datasets}
