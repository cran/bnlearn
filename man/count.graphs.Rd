\name{graph enumeration}
\alias{count.graphs}
\alias{graph enumeration}
\title{Count graphs with specific characteristics}
\description{
  Count directed acyclic graphs of various sizes and/or with specific
  characteristics.
}
\usage{
count.graphs(type = "all.dags", ..., debug = FALSE)
}
\arguments{
  \item{type}{a character string, the label describing the types of graphs to
    be counted (see below).}
  \item{\dots}{additional parameters (see below).}
  \item{debug}{a boolean value. If \code{TRUE}, a lot of debugging output is
    printed. Otherwise, the function is completely silent. Ignored in some
    generation methods.}
}
\details{

  The types of graphs, and the associated additional parameters, are:

  \itemize{

    \item \code{"all-dags"}: all directed acyclic graphs with \code{nodes}
      nodes.
    \item \code{"dags-given-ordering"}: all directed acyclic graphs
      with a specific topological ordering and \code{nodes} nodes.
    \item \code{"dags-with-k-roots"}: all directed acyclic graphs with \code{k}
      root nodes and \code{nodes} nodes.
    \item \code{"dags-with-r-arcs"}: all directed acyclic graphs with \code{r}
      arcs and \code{nodes} nodes.
    \item \code{"dags-in-equivalence-class"}: all directed acyclic arcs in the
      equivalence class described by \code{eqclass}.

  }

}
\value{

  \code{count.graphs()} returns an object of class \code{bigz} from the
  \pkg{gmp} package, a vector with the graph counts.

}
\references{

  Harary F, Palmer EM (1973). "Graphical Enumeration." Academic Press.

  Rodionov VI (1992). "On the Number of Labeled Acyclic Digraphs."
    \emph{Discrete Mathematics}, \strong{105}:319--321.

  Liskovets VA (1976). "On the Number of Maximal Vertices of a Random Acyclic
    Digraph." \emph{Theory of Probability and its Applications},
    \strong{20}(2):401--409.

  Wienobst M, Luttermann M, Bannach M, Liskiewicz (2023). "Efficient Enumeration
    of Markov Equivalent DAGs." \emph{Proceedings of the 37th AAAI Conference on
    Artificial Intelligence}, 12313--12320.

}
\examples{
\dontrun{
count.graphs("dags.with.r.arcs", nodes = 3:6, r = 2)
}}
\author{Marco Scutari}
\keyword{graphs}
