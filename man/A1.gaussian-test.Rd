\name{gaussian.test}
\docType{data}
\alias{gaussian.test}
\title{ Synthetic (continuous) dataset to test learning algorithms }
\description{

  This a synthetic dataset used as a test case in the
  \pkg{bnlearn} package.

}
\usage{
  gaussian.test
}
\details{

  The \code{gaussian.test} dataset containts the seven normal
  (gaussian) variables.

}
\examples{

# these are the R commands used to generate this dataset.
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
