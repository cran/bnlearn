\name{compare}
\alias{compare}
\alias{all.equal.bn}
\alias{shd}
\alias{hamming}
\title{Compare two different Bayesian networks}
\description{

  Compare two different Bayesian networks; compute their Structural Hamming
  Distance (SHD) or the Hamming distance between their skeletons.

}
\usage{
compare(target, current, arcs = FALSE)
\method{all.equal}{bn}(target, current, ...)

shd(learned, true, wlbl = FALSE, debug = FALSE)
hamming(learned, true, debug = FALSE)
}
\arguments{
  \item{target, learned}{an object of class \code{bn}.}
  \item{current, true}{another object of class \code{bn}.}
  \item{\dots}{extra arguments from the generic method (currently ignored).}
  \item{wlbl}{a boolean value. If \code{TRUE} arcs whose directions have been
    fixed by a whitelist or a by blacklist are preserved when constructing
    the CPDAGs of \code{learned} and \code{true}.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
  \item{arcs}{a boolean value. See below.}
}
\note{

  Note that SHD, as defined in the reference, is defined on CPDAGs; therefore
  \code{cpdag} is called on both \code{learned} and \code{true} before computing
  the distance.

}
\value{

  \code{compare} returns a list containing the number of true positives
  (\code{tp}, the number of arcs in \code{current} also present in
  \code{target}), of false positives (\code{fp}, the number of arcs in
  \code{current} not present in \code{target}) and of false negatives
  (\code{fn}, the number of arcs not in \code{current} but present in
  \code{target}) if \code{arcs} is \code{FALSE}; or the corresponding arc sets
  if \code{arcs} is \code{TRUE}.

  \code{all.equal} returns either \code{TRUE} or a character string describing
  the differences between \code{target} and \code{current}.

  \code{shd} and \code{hamming} return a non-negative integer number.

}
\examples{
data(learning.test)

e1 = model2network("[A][B][C|A:B][D|B][E|C][F|A:E]")
e2 = model2network("[A][B][C|A:B][D|B][E|C:F][F|A]")
shd(e2, e1, debug = TRUE)
unlist(compare(e1,e2))
compare(target = e1, current = e2, arcs = TRUE)
}
\references{

  Tsamardinos I, Brown LE, Aliferis CF (2006). "The Max-Min Hill-Climbing
      Bayesian Network Structure Learning Algorithm". \emph{Machine Learning},
      \strong{65}(1), 31-78.

}
\author{ Marco Scutari }
\keyword{graphs}
