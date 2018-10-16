\name{graph enumeration}
\alias{count.graphs}
\alias{graph enumeration}
\title{Count graphs with specific characteristics}
\description{
  Count directed acyclic graphs of various sizes with specific characteristics.
}
\usage{
count.graphs(type = "all.dags", nodes, ..., debug = FALSE)
}
\arguments{
  \item{type}{a character string, the label describing the types of graphs to
    be counted (see below).}
  \item{nodes}{a vector of positive integers, the graph sizes as given by the
    numbers of nodes.}
  \item{\dots}{additional parameters (see below).}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent. Ignored in some
    generation methods.}
}
\details{

  The types of graphs, and the associated additional parameters, are:

  \itemize{

    \item \code{all-dags}: all directed acyclic graphs.
    \item \code{dags-given-ordering}: all directed acyclic graphs
      with a specific topological ordering.
    \item \code{dags-with-k-roots}: all directed acyclic graphs with \code{k}
      root nodes.
    \item \code{dags-with-r-arcs}: all directed acyclic graphs with \code{r}
      arcs.

  }

}
\value{

  \code{count.graphs()} returns an objects of class \code{bigz} from the
  \pkg{gmp} package, a vector with the graph counts.

}
\references{

  Harary F, Palmer EM (1973). "Graphical Enumeration". Academic Press.

  Rodionov VI (1992). "On the Number of Labeled Acyclic Digraphs".
    \emph{Discrete Mathematics}, \strong{105}:319--321.

  Liskovets VA (1976). "On the Number of Maximal Vertices of a Random Acyclic
    Digraph". \emph{Theory of Probability and its Applications},
    \strong{20}(2):401--409.

}
\examples{
\dontrun{
count.graphs("dags.with.r.arcs", nodes = 3:6, r = 2)
}}
\author{Marco Scutari}
\keyword{graphs}
