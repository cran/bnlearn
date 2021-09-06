\name{strength.plot}
\alias{strength.plot}
\title{Arc strength plot}
\description{

  Plot a Bayesian network and format its arcs according to the strength of the
  dependencies they represent. Requires the \pkg{Rgraphviz} package.

}
\usage{
strength.plot(x, strength, threshold, cutpoints, highlight = NULL, groups,
  layout = "dot", shape = "circle", main = NULL, sub = NULL, render = TRUE,
  debug = FALSE)
}
\arguments{
  \item{x}{an object of class \code{bn}.}
  \item{strength}{an object of class \code{bn.strength} computed from the object
    of class \code{bn} corresponding to the \code{x} argument.}
  \item{threshold}{a numeric value. See below.}
  \item{cutpoints}{an array of numeric values. See below.}
  \item{highlight}{a list, see \code{\link{graphviz.plot}} for details.}
  \item{groups}{a list of character vectors, representing groups of node labels
    of nodes that should be plotted close to each other.}
  \item{layout}{a character string, the layout argument that will be passed to
    \pkg{Rgraphviz}. Possible values are \code{dots}, \code{neato},
    \code{twopi}, \code{circo} and \code{fdp}. See \pkg{Rgraphviz} documentation
    for details.}
  \item{shape}{a character string, the shape of the nodes. Can be \code{circle},
    \code{ellipse} or \code{rectangle}.}
  \item{main}{a character string, the main title of the graph. It's plotted at
    the top of the graph.}
  \item{sub}{a character string, a subtitle which is plotted at the bottom of
    the graph.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
  \item{render}{a logical value. If \code{TRUE}, \code{strength.plot()} actually
    draws the figure in addition to returning the corresponding \code{graph}
    object. If \code{FALSE}, no figure is produced.}
}
\details{

  The \code{threshold} argument is used to determine which arcs are supported
  strongly enough by the data to be deemed significant:
  \itemize{
    \item if arc strengths have been computed using conditional independence
      tests, any strength coefficient (which is the p-value of the test) lesser
      than or equal to the threshold is considered significant.
    \item if arc strengths have been computed using network scores, any strength
      coefficient (which is the increase/decrease of the network score caused by
      the removal of the arc) lesser than the threshold is considered
      significant.
    \item if arc strengths have been computed using bootstrap or using Bayes
      factors, any strength coefficient (which can be interpreted as a
      probability for inclusion) greater or equal than the threshold is
      considered significant.
  }
  The default value is the value of the \code{strength} attribute of the
  \code{bn.strength} object passed via the \code{strength} argument.

  Non-significant arcs are plotted as dashed lines.

  The \code{cutpoints} argument is an array of numeric values used to divide
  the range of the strength coefficients into intervals. The interval each
  strength coefficient falls into determines the line width of the corresponding
  arc in the plot. The default intervals are delimited by

    \code{  unique(c(0, threshold/c(10, 5, 2, 1.5, 1), 1))}

  if the coefficients are computed from conditional independence tests, by

    \code{  unique(c(0, (1 - threshold)/c(10, 5, 2, 1.5, 1), 1))}

  for bootstrap estimates or by the quantiles

    \code{  quantile(s[s < threshold], 1 - c(0.50, 0.75, 0.90, 0.95, 1))}

  of the significant differences if network scores are used.

}
\value{

  \code{graphviz.plot()} returns invisibly the \code{graph} object produced by
  \pkg{Rgraphviz}. It can be further modified using the commands present in the
  \pkg{graph} and \pkg{Rgraphviz} packages, and it contains the arc strengths
  in the edge weight attribute.

}
\examples{
\dontrun{
# plot the network learned by hc().
res = hc(learning.test)
strength = arc.strength(res, learning.test, criterion = "x2")
\dontrun{strength.plot(res, strength)}
# add another (non-significant) arc and plot the network again.
res = set.arc(res, "A", "C")
strength = arc.strength(res, learning.test, criterion = "x2")
\dontrun{strength.plot(res, strength)}
}
}
\author{Marco Scutari}
\keyword{plots}
\keyword{interfaces to other packages}
