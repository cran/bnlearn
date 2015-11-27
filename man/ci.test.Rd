\name{ci.test}
\alias{ci.test}
\title{Independence and Conditional Independence Tests}
\description{

  Perform either an independence test or a conditional independence test.

}
\usage{
ci.test(x, y, z, data, test, B, debug = FALSE)
}
\arguments{
  \item{x}{a character string (the name of a variable), a data frame, a numeric
    vector or a factor object.}
  \item{y}{a character string (the name of another variable), a numeric vector
    or a factor object.}
  \item{z}{a vector of character strings (the names of the conditioning
    variables), a numeric vector, a factor object or a data frame. If
    \code{NULL} an independence test will be executed.}
  \item{data}{a data frame containing the variables to be tested.}
  \item{test}{a character string, the label of the conditional independence
    test to be used in the algorithm. If none is specified, the default test
    statistic is the \emph{mutual information} for categorical variables, the
    Jonckheere-Terpstra test for ordered factors and the \emph{linear 
    correlation} for continuous variables. See \code{\link{bnlearn-package}}
    for details.}
  \item{B}{a positive integer, the number of permutations considered for each
    permutation test. It will be ignored with a warning if the conditional
    independence test specified by the \code{test} argument is not a 
    permutation test.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
}
\value{

  An object of class \code{htest} containing the following components:

  \item{statistic}{the value the test statistic.}
  \item{parameter}{the degrees of freedom of the approximate chi-squared or t
    distribution of the test statistic; the number of permutations computed by
    Monte Carlo tests. Semiparametric tests have both.}
  \item{p.value}{the p-value for the test.}
  \item{method}{a character string indicating the type of test performed, and
    whether Monte Carlo simulation or continuity correction was used.}
  \item{data.name}{a character string giving the name(s) of the data.}
  \item{null.value}{the value of the test statistic under the null hypothesis,
    always 0.}
  \item{alternative}{a character string describing the alternative hypothesis.}

}
\references{

  \bold{for parametric and discrete permutation tests:}

  Edwards DI (2000). \emph{Introduction to Graphical Modelling}.
      Springer, 2nd edition.

  \bold{for shrinkage tests:}

  Hausser J, Strimmer K (2009). "Entropy inference and the James-Stein
      estimator, with application to nonlinear gene association networks".
      \emph{Statistical Applications in Genetics and Molecular Biology},
      \strong{10}, 1469-1484.

  Ledoit O, Wolf M (2003). "Improved Estimation of the Covariance Matrix
      of Stock Returns with an Application to Portfolio Selection".
      \emph{Journal of Empirical Finance}, \strong{10}, 603-621.

  \bold{for continuous permutation tests:}

  Legendre P (2000). "Comparison of Permutation Methods for the Partial
      Correlation and Partial Mantel Tests". \emph{Journal of Statistical
      Computation and Simulation}, \strong{67}, 37-73.

  \bold{for semiparametric discrete tests:}

  Tsamardinos I, Borboudakis G (2010). "Permutation Testing Improves Bayesian
      Network Learning". In "Machine Learning and Knowledge Discovery in
      Databases", pp. 322-337. Springer.

}
\examples{
data(gaussian.test)
data(learning.test)

# using a data frame and column labels.
ci.test(x = "F" , y = "B", z = c("C", "D"), data = gaussian.test)
# using a data frame.
ci.test(gaussian.test)
# using factor objects.
attach(learning.test)
ci.test(x = F , y = B, z = data.frame(C, D))
}
\seealso{\code{\link{choose.direction}}, \code{\link{arc.strength}}.}
\author{ Marco Scutari }
\keyword{htest}
