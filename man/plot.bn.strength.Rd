\name{plot.bn.strength}
\alias{plot.bn.strength}
\title{Plot arc strengths derived from bootstrap}
\description{

  Plot arc strengths derived from bootstrap resampling.

}
\usage{
\method{plot}{bn.strength}(x, draw.threshold = TRUE, main = NULL,
  xlab = "arc strengths", ylab = "CDF(arc strengths)", ...)
}
\arguments{
  \item{x}{an object of class \code{bn.strength}.}
  \item{draw.threshold}{a boolean value. If \code{TRUE}, a dashed vertical line
    is drawn at the threshold.}
  \item{main,xlab,ylab}{character strings, the main title and the axes labels.}
  \item{\dots}{other graphical parameters.}
}
\note{
  The \code{xlim} and \code{ylim} graphical parameters are always overridden.
}
\examples{
data(learning.test)

start = random.graph(nodes = names(learning.test), num = 50)
netlist = lapply(start, function(net) {
  hc(learning.test, score = "bde", iss = 10, start = net) })
arcs = custom.strength(netlist, nodes = names(learning.test), cpdag = FALSE)
plot(arcs)
}
\author{ Marco Scutari }
\keyword{hplot}
