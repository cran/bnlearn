\name{ci.test}
\alias{ci.test}
\title{Independence and conditional independence tests}
\description{

  Perform an independence or a conditional independence test.

}
\usage{
ci.test(x, y, z, data, test, ..., debug = FALSE)
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
    correlation} for continuous variables. See \code{\link{independence tests}}
    for details.}
  \item{\dots}{optional arguments to be passed to the test specified by
    \code{test}. See below for details.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
}
\details{

  Additional arguments of the \code{ci.test()} function:

  \itemize{

    \item \code{B}: a positive integer, the number of permutations used to
      compute the p-value of permutation tests. The default value is \code{5000}
      for nonparametric permutation tests and \code{100} for semiparametric
      permutation tests.
    \item \code{fun}: the function that computes the conditional independence
      test in the \code{custom-test} test. \code{fun} must have arguments
      \code{x}, \code{y}, \code{z}, \code{data} and \code{args}, in this order;
      in other words, it must have signature \code{function(x, y, z, data,
      args)}. \code{x} and \code{y} will be the labels of the nodes to test for
      independence (one character string each); \code{z} will be the labels of
      the nodes in the conditioning set (a vector of character strings,
      possibily \code{NULL} for empty sets); \code{data} will contain the
      complete data set, with all the variables (a data frame); and \code{args}
      will be a list containing any additional arguments to the test.
    \item \code{args}: a list containing the optional arguments to \code{fun},
      for tuning \code{custom-test} test functions.

  }

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
\seealso{\code{\link{independence tests}}, \code{\link{arc.strength}}.}
\author{Marco Scutari}
\keyword{independence tests}
