\name{graphviz.plot}
\alias{graphviz.plot}
\title{ Advanced Bayesian network plots }
\description{

  Plot the graph associated with a Bayesian network
  using the \pkg{Rgraphviz} package.

}
\usage{
  graphviz.plot(x, highlight = NULL, layout = "dot")
}
\arguments{
  \item{x}{an object of class \code{bn}.}
  \item{highlight}{a list, see below.}
  \item{layout}{the layout parameter to be passed to
    \pkg{Rgraphviz}. Possible values are \code{dots}, 
    \code{neato}, \code{twopi}, \code{circo} and 
    \code{fdp}.}
}
\details{

  The highlight parameter is a list with at least
  one of the following elements:

  \itemize{

    \item \code{nodes}: a character vector, the labels
      of the nodes to be highlighted.
    \item \code{arcs}: the arcs to be highlighted (a 
      two-column matrix, whose columns are labeled 
      \code{from} and \code{to}).

  }

  and optionally one or more of the following formatting 
  parameters:

  \itemize{

    \item \code{col}: an integer or character string (the 
      highlight colour). The default value is \code{red}.

    \item \code{fill}: an integer or character string (the 
      colour used as a background colour for the nodes).
      The default value is \code{white}.

  }

}
\author{ Marco Scutari }
\keyword{hplot}
