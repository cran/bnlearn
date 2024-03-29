\name{graphviz.plot}
\alias{graphviz.plot}
\title{Advanced Bayesian network plots}
\description{

  Plot the graph associated with a Bayesian network using the \pkg{Rgraphviz}
  package.

}
\usage{
graphviz.plot(x, highlight = NULL, groups, layout = "dot",
  shape = "rectangle", fontsize = 12, main = NULL, sub = NULL,
  render = TRUE)
}
\arguments{
  \item{x}{an object of class \code{bn} or \code{bn.fit}.}
  \item{highlight}{a list, see below.}
  \item{groups}{a list of character vectors, representing groups of node labels
    of nodes that should be plotted close to each other.}
  \item{layout}{a character string, the layout argument that will be passed to
    \pkg{Rgraphviz}. Possible values are \code{dots}, \code{neato},
    \code{twopi}, \code{circo} and \code{fdp}. See \pkg{Rgraphviz}
    documentation for details.}
  \item{shape}{a character string, the shape of the nodes. Can be \code{circle},
    \code{ellipse} or \code{rectangle}.}
  \item{fontsize}{a positive number, the font size for the node labels.}
  \item{main}{a character string, the main title of the graph. It's plotted at
    the top of the graph.}
  \item{sub}{a character string, a subtitle which is plotted at the bottom of
    the graph.}
  \item{render}{a logical value. If \code{TRUE}, \code{graphviz.plot()} actually
    draws the figure in addition to returning the corresponding \code{graph}
    object. If \code{FALSE}, no figure is produced.}
}
\details{

  The \code{highlight} argument is a list with at least one of the following
  elements:

  \itemize{

    \item \code{nodes}: a character vector, the labels of the nodes to be
      highlighted.
    \item \code{arcs}: the arcs to be highlighted (a two-column matrix, whose
      columns are labeled \code{from} and \code{to}).

  }

  and optionally one or more of the following graphical parameters:

  \itemize{

    \item \code{col}: an integer or character string (the highlight colour for
      the arcs and the node frames). The default value is \code{red}.
    \item \code{textCol}: an integer or character string (the highlight colour
      for the labels of the nodes). The default value is \code{black}.
    \item \code{fill}: an integer or character string (the colour used as a
      background colour for the nodes). The default value is \code{transparent}.
    \item \code{lwd}: a positive number (the line width of highlighted arcs).
      It overrides the line width settings in \code{strength.plot()}. The
      default value is to use the global settings of \pkg{Rgraphviz}.
    \item \code{lty}: the line type of highlighted arcs. Possible values are
      0, 1, 2, 3, 4, 5, 6, "blank", "solid", "dashed", "dotted", "dotdash",
      "longdash" and "twodash". The default value is to use the global
      settings of \pkg{Rgraphviz}.

  }

  Note that all these parameters take a single value that is then applied to all
  nodes and arcs that will be highlighted.

}
\value{

  \code{graphviz.plot()} returns invisibly the \code{graph} object produced from
  the network passed as the \code{x} argument. It can be further modified using
  the \pkg{graph} and \pkg{Rgraphviz} packages.

}
\seealso{\code{\link{plot.bn}}.}
\author{Marco Scutari}
\keyword{interfaces to other packages}
\keyword{plots}
