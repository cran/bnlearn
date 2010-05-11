\name{strength.plot}
\alias{strength.plot}
\title{ Arc strength plot }
\description{

  Plot a Bayesian network and format its arcs according to the
  strength of the dependencies they represent. Requires the
  \pkg{Rgraphviz} package.

}
\usage{
strength.plot(x, strength, threshold, cutpoints, highlight = NULL,
  layout = "dot", shape = "circle", main = NULL, sub = NULL,
  debug = FALSE)
}
\arguments{
  \item{x}{an object of class \code{bn}.}
  \item{strength}{an object of class \code{bn.strength} computed from
    the object of class \code{bn} corresponding to the \code{x}
    parameter.}
  \item{threshold}{a numeric value. See below.}
  \item{cutpoints}{an array of numeric values. See below.}
  \item{highlight}{a list, see \code{\link{graphviz.plot}} for details.}
  \item{layout}{a character string, the layout parameter to
    be passed to \pkg{Rgraphviz}. Possible values are \code{dots},
    \code{neato}, \code{twopi}, \code{circo} and \code{fdp}. See
    \pkg{Rgraphviz} documentation for details.}
  \item{shape}{a character string, the shape of the nodes. Can be
    either \code{circle} or \code{ellipse}.}
  \item{main}{a character string, the main title of the graph.
    It's plotted at the top of the graph.}
  \item{sub}{a character string, a subtitle which is plotted
    at the bottom of the graph.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output
       is printed; otherwise the function is completely silent.}
}
\details{

  The \code{threshold} parameter is used to determine which arcs are supported
  strongly enough by the data to be deemed significant:
  \itemize{
    \item if arc strengths have been computed using conditional independence tests,
      any strength coefficient (which is the p-value of the test) lesser
      or equal than the threshold is considered significant. In this case
      the default value of \code{threshold} is equal to the value of the
      \code{alpha} parameter used in the call to \code{arc.strength},
      which in turn defaults to the one used by the learning algorithm
      (if any) or to \code{0.05}.
    \item if arc strengths have been computed using network scores,
      any strength coefficient (which is the increase/decrease of the network
      score caused by the removal of the arc) lesser than the threshold is
      considered significant. In this case the default value of
      \code{threshold} is \code{0}.
    \item if arc strengths have been computed using bootstrap, any
      strength coefficient (which is the relative frequency of the
      arc in the networks learned from the bootstrap replicates) greater
      or equal than the threshold is considered significant. In this case
      the default value of \code{threshold} is \code{0.5}.
  }

  Non-significant arcs are plotted as dashed lines.

  The \code{cutpoints} parameter is an array of numeric values used to
  divide the range of the strength coefficients into intervals. The interval
  each strength coefficient falls into determines the line width of the
  corresponding arc in the plot. The default intervals are delimited by

    \code{  unique(c(0, threshold/c(10, 5, 2, 1.5, 1), 1))}

  if the coefficients are computed from conditional independence tests, by

    \code{  1 - unique(c(0, threshold/c(10, 5, 2, 1.5, 1), 1))}

  for bootstrap estimates or by the quantiles

    \code{  quantile(-s[s < threshold], c(0.50, 0.75, 0.90, 0.95, 1))}

  of the significant coefficients if network scores are used.

}
\value{

  The object of class \code{graphAM} used to format and render the
  plot. It can be further modified using the commands present in the
  \pkg{graph} and \pkg{Rgraphviz} packages.

}
\examples{
\dontrun{
# plot the network learned by gs().
res = set.arc(gs(learning.test), "A", "B")
strength = arc.strength(res, learning.test, criterion = "x2")
strength.plot(res, strength)
# add another (non-significant) arc and plot the network again.
res = set.arc(res, "A", "C")
strength = arc.strength(res, learning.test, criterion = "x2")
strength.plot(res, strength)
}
}
\author{ Marco Scutari }
\keyword{hplot}
