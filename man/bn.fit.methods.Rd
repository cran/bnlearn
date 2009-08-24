\name{bn.fit utilities}
\alias{bn.fit utilities}
\alias{fitted.bn.fit}
\alias{coef.bn.fit}
\alias{residuals.bn.fit}
\alias{fitted.bn.fit.gnode}
\alias{coef.bn.fit.gnode}
\alias{residuals.bn.fit.gnode}
\alias{fitted.bn.fit.dnode}
\alias{coef.bn.fit.dnode}
\alias{residuals.bn.fit.dnode}

\title{ Utilities to manipulate fitted Bayesian networks }
\description{

 Assign or extract various quantities of interest from an
 object of class \code{bn.fit}, \code{bn.fit.dnode} or
 \code{bn.fit.gnode}.

}
\usage{
  ## methods available for "bn.fit"
  \method{fitted}{bn.fit}(object, ...)
  \method{coef}{bn.fit}(object, ...)
  \method{residuals}{bn.fit}(object, ...)

  ## methods available for "bn.fit.dnode"
  \method{coef}{bn.fit.gnode}(object, ...)

  ## methods available for "bn.fit.gnode"
  \method{fitted}{bn.fit.gnode}(object, ...)
  \method{coef}{bn.fit.gnode}(object, ...)
  \method{residuals}{bn.fit.gnode}(object, ...)
}
\arguments{
  \item{object}{an object of class \code{bn.fit}, \code{bn.fit.dnode}
    or \code{bn.fit.gnode}.}
  \item{\dots}{additional arguments (currently ignored).}
}
\details{

  \code{coef} (and its alias \code{coefficients}) extracts model
  coefficients (which are conditional probabilities in discrete
  networks and linear regression coefficients in Gaussian networks).

  \code{residuals} (and its alias \code{resid}) extracts model
  residuals and \code{fitted} (and its alias \code{fitted.values})
  extracts fitted values from fitted Gaussian networks.

}
\value{

  A list with an element for each node in the network (if \code{object}
  has class \code{bn.fit}) or a numeric vector (if \code{object} has class
  \code{bn.fit.dnode} or \code{bn.fit.gnode}).

}
\seealso{
  \code{\link{bn.fit}}, \code{\link{bn.fit-class}}.
}
\examples{
data(gaussian.test)
res = hc(gaussian.test)
fitted = bn.fit(res, gaussian.test)

coefficients(fitted)
# $A
# (Intercept)
#    1.007493
#
# $B
# (Intercept)
#    2.039499
#
# $C
# (Intercept)           A           B
#    2.001083    1.995901    1.999108
#
# $D
# (Intercept)           B
#    5.995036    1.498395
#
# $E
# (Intercept)
#    3.493906
#
# $F
#  (Intercept)            A            D            E            G
# -0.006047321  1.994853041  1.005636909  1.002577002  1.494373265
#
# $G
# (Intercept)
#    5.028076
#
coefficients(fitted$C)
# (Intercept)           A           B
#    2.001083    1.995901    1.999108
str(residuals(fitted))
# List of 7
#  $ A: num [1:5000] 0.106 -1.255 0.847 -0.174 -0.519 ...
#  $ B: num [1:5000] -0.107 9.295 0.993 1.818 2.473 ...
#  $ C: num [1:5000] -1.01 0.183 -0.677 -0.153 -1.997 ...
#  $ D: num [1:5000] -0.23 0.377 0.518 0.162 -0.22 ...
#  $ E: num [1:5000] -2.612 3.546 0.341 -2.488 0.591 ...
#  $ F: num [1:5000] -0.861 1.271 -0.262 -0.479 -0.782 ...
#  $ G: num [1:5000] 4.1883 -1.3492 -2.6036 1.0574 0.0895 ...

data(learning.test)
res2 = hc(learning.test)
fitted2 = bn.fit(res2, learning.test)

coefficients(fitted2$E)
# , , F = a
#
#    B
# E        a      b      c
#   a 0.1902 0.0126 0.0244
#   b 0.0230 0.0110 0.0234
#   c 0.0230 0.0376 0.1566
#
# , , F = b
#
#    B
# E        a      b      c
#   a 0.0946 0.0166 0.0498
#   b 0.1158 0.0192 0.1062
#   c 0.0258 0.0166 0.0536
}
\author{ Marco Scutari }
\keyword{utilities}
\keyword{graphs}
