\name{arc operations}
\alias{arc operations}
\alias{set.arc}
\alias{drop.arc}
\alias{reverse.arc}
\alias{set.edge}
\alias{drop.edge}
\title{Drop, add or set the direction of an arc or an edge}
\description{

  Drop, add or set the direction of a directed or undirected arc (also
  known as edge).

}
\usage{
# arc operations.
set.arc(x, from, to, check.cycles = TRUE, check.illegal = TRUE, debug = FALSE)
drop.arc(x, from, to, debug = FALSE)
reverse.arc(x, from, to, check.cycles = TRUE, check.illegal = TRUE, debug = FALSE)

# edge (i.e. undirected arc) operations
set.edge(x, from, to, check.cycles = TRUE, check.illegal = TRUE, debug = FALSE)
drop.edge(x, from, to, debug = FALSE)
}
\arguments{
  \item{x}{an object of class \code{bn}.}
  \item{from}{a character string, the label of a node.}
  \item{to}{a character string, the label of another node.}
  \item{check.cycles}{a boolean value. If \code{TRUE} the graph is tested for
    acyclicity; otherwise the graph is returned anyway.}
  \item{check.illegal}{a boolean value. If \code{TRUE} arcs that break the
    parametric assumptions of \code{x}, such as those from continuous to
    discrete nodes in conditional Gaussian networks, cause an error.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
}
\details{

  The \code{set.arc()} function operates in the following way:
  \itemize{
    \item if there is no arc between \code{from} and \code{to}, the arc
      \code{from} \eqn{\rightarrow}{->} \code{to} is added.
    \item if there is an undirected arc between \code{from} and \code{to}, its
      direction is set to \code{from} \eqn{\rightarrow}{->} \code{to}.
    \item if the arc \code{to} \eqn{\rightarrow}{->} \code{from} is present,
      it is reversed.
    \item if the arc \code{from} \eqn{\rightarrow}{->} \code{to} is present, no
      action is taken.
  }

  The \code{drop.arc()} function operates in the following way:
  \itemize{
    \item if there is no arc between \code{from} and \code{to}, no action is
      taken.
    \item if there is a directed or an undirected arc between \code{from} and
      \code{to}, it is dropped regardless of its direction.
  }

  The \code{reverse.arc()} function operates in the following way:
  \itemize{
    \item if there is no arc between \code{from} and \code{to}, it returns an
      error.
    \item if there is an undirected arc between \code{from} and \code{to}, it
      returns an error.
    \item if the arc \code{to} \eqn{\rightarrow}{->} \code{from} is present, it
      is reversed.
    \item if the arc \code{from} \eqn{\rightarrow}{->} \code{to} is present, it
      is reversed.
  }

  The \code{set.edge()} function operates in the following way:
  \itemize{
    \item if there is no arc between \code{from} and \code{to}, the undirected
      arc \code{from} - \code{to} is added.
    \item if there is an undirected arc between \code{from} and \code{to}, no
      action is taken.
    \item if either the arc \code{from} \eqn{\rightarrow}{->} \code{to} or the
      arc \code{to} \eqn{\rightarrow}{->} \code{from} are present, they are
      replaced with the undirected arc \code{from} - \code{to}.
  }

  The \code{drop.edge()} function operates in the following way:
  \itemize{
    \item if there is no undirected arc between \code{from} and \code{to}, no
      action is taken.
    \item if there is an undirected arc between \code{from} and \code{to}, it
      is removed.
    \item if there is a directed arc between \code{from} and \code{to}, no
      action is taken.
  }

}
\value{

  All functions return invisibly an updated copy of \code{x}.

}
\examples{
dag = cpdag(model2network("[A][C][F][B|A][D|A:C][E|B:F]"))
dag

## use debug = TRUE to get more information.
updated = set.arc(dag, "A", "B")
updated
updated = drop.arc(dag, "A", "B")
updated
updated = drop.edge(dag, "A", "B")
updated
updated = reverse.arc(dag, "A", "D")
updated
}
\author{Marco Scutari}
\keyword{graphs}
