\name{ci.test}
\alias{ci.test}
\alias{ci.test.character}
\alias{ci.test.data.frame}
\alias{ci.test.numeric}
\alias{ci.test.factor}
\alias{ci.test.default}
\title{ Independence and Conditional Independence Tests }
\description{

  Perform an independence or a  conditional independence test.

}
\usage{
  \method{ci.test}{character}(x, y = NULL, z = NULL, data, test = NULL,
    debug = FALSE, ...)
  \method{ci.test}{data.frame}(x, test = NULL, debug = FALSE, ...)
  \method{ci.test}{numeric}(x, y = NULL, z = NULL, test = NULL,
    debug = FALSE, ...)
  \method{ci.test}{factor}(x, y = NULL, z = NULL, test = NULL,
    debug = FALSE, ...)
  \method{ci.test}{default}(x, ...)
}
\arguments{
  \item{x}{a character string (the name of a variable), a data frame, a numeric
      vector or a factor object.}
  \item{y}{a character string (the name of another variable), a numeric
      vector or a factor object.}
  \item{z}{a vector of character strings (the names of the conditioning
      variables), a numeric vector, a factor object or a data frame.
      If \code{NULL} an independence test will be executed.}
  \item{data}{a data frame, containing the variables to be tested.}
  \item{test}{a character string, the label of the conditional
      independence test to be used in the algorithm. Possible
      values are \code{mi} (\emph{mutual information} for discrete
      data), \code{fmi} (\emph{fast mutual information}), \code{aict}
      (\emph{AIC}), \code{x2} (\emph{Pearson's \eqn{\chi^2}{chi-square}}),
      \code{cor} (\emph{linear correlation}), \code{zf} (\emph{Fisher's Z})
      and \code{mi-g} (\emph{mutual information} for gaussian data).
      If none is specified, the default one is \emph{mutual
      information} for discrete data sets and \emph{linear correlation}
      for continuous ones. See \code{\link{bnlearn-package}} for details.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output
       is printed; otherwise the function is completely silent.}
  \item{\dots}{extra arguments from the generic method (currently ignored).}
}
\value{

  An object of class \code{htest} containing the following components:

  \item{statistic}{the value the conditional independence test statistic.}
  \item{parameter}{the degrees of freedom of the approximate
    chi-squared or t distribution of the test statistic, \code{NA} if the
    p-value is computed by Monte Carlo simulation.}
  \item{p.value}{the p-value for the test.}
  \item{method}{a character string indicating the type of test
    performed, and whether Monte Carlo simulation or continuity
    correction was used.}
  \item{data.name}{a character string giving the name(s) of the data.}
  \item{null.value}{the value of the test statistic under the null
    hypothesis, always 0.}
  \item{alternative}{a character string describing the alternative
    hypothesis}

}
\examples{
data(gaussian.test)
data(learning.test)

# using a data frame and column labels.
ci.test(x = "F" , y = "B", z = c("C", "D"), data = gaussian.test)
#
#	 linear correlation
#
# data:  F ~ B | C + D
# cor = -0.0423, df = 4996, p-value = 0.002797
# alternative hypothesis: true value is not equal to 0

# using a data frame.
ci.test(gaussian.test)
#
#	 linear correlation
#
# data:  A ~ B | C + D + E + F + G
# cor = -0.1535, df = 4993, p-value < 2.2e-16
# alternative hypothesis: true value is not equal to 0

# using factor objects.
ci.test(x = learning.test[, "F"] , y = learning.test[, "B"],
  z = learning.test[, c("C", "D")] )
#
#	 mutual information
#
# data:  learning.test[, "F"] ~ learning.test[, "B"] | learning.test[, c("C", "D")]
# mi = 25.2664, df = 18, p-value = 0.1178
# alternative hypothesis: true value is greater than 0

}
\seealso{\code{\link{choose.direction}}, \code{\link{arc.strength}}. }
\author{ Marco Scutari }
\keyword{htest}
