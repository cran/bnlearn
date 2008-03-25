\name{arc operations}
\alias{set.arc}
\alias{drop.arc}
\alias{reverse.arc}
\title{ Drop, add or set the direction of an arc }
\description{

  Drop, add or set the direction of an arc.

}
\usage{
  set.arc(x, from, to, check.cycles = TRUE, debug = FALSE)
  drop.arc(x, from, to, debug = FALSE)
  reverse.arc(x, from, to, check.cycles = TRUE, debug = FALSE)
}
\arguments{
  \item{x}{an object of class "bn".}
  \item{from}{a character string, the label of a node.}
  \item{to}{a character string, the label of another node.}
  \item{check.cycles}{a boolean value. If \code{TRUE} the graph is
    tested for acyclicity; otherwise the graph is returned anyway.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output
       is printed; otherwise the function is completely silent.}
}
\details{

  The \code{set.arc} function operates in the following way:
  \itemize{
    \item if there is no arc between \code{from} and \code{to},
      the arc (\code{from}, \code{to}) is added.
    \item if there is an undirected arc between \code{from} and
      \code{to}, its direction is set to (\code{from}, \code{to}).
    \item if the arc (\code{to}, \code{from}) is present, it's
      reversed.
    \item if the arc (\code{from}, \code{to}) is present, no
      action is taken.
  }

  The \code{drop.arc} function operates in the following way:
  \itemize{
    \item if there is no arc between \code{from} and \code{to},
      no action is taken.
    \item if there is an undirected arc between \code{from} and
      \code{to}, it's dropped.
    \item if there is a directed arc between \code{from} and
      \code{to}, it's dropped regardless of its direction.
  }

  The \code{reverse.arc} function operates in the following way:
  \itemize{
    \item if there is no arc between \code{from} and \code{to},
      it returns an error.
    \item if there is an undirected arc between \code{from} and
      \code{to}, it returns an error.
    \item if the arc (\code{to}, \code{from}) is present, it's
      reversed.
    \item if the arc (\code{from}, \code{to}) is present, it's
      reversed.
  }

}
\value{

  \code{set.arc} and \code{drop.arc} return invisibly an updated
  copy of \code{x}.

}
\examples{
data(learning.test)

res = gs(learning.test)
res2 = set.arc(res, "A", "B")
compare(res, res2, debug = TRUE)
# * children of A in r1 not present in r2:
# character(0)
# * children of A in r2 not present in r1:
# [1] "B"
# * parents of B in r1 not present in r2:
# character(0)
# * parents of B in r2 not present in r1:
# [1] "A"
# * directed arcs in r1 not present in r2:
# character(0)
# * directed arcs in r2 not present in r1:
# [1] "A -> B"
# * undirected arcs in r1 not present in r2:
# [1] "A - B" "B - A"
# * undirected arcs in r2 not present in r1:
# character(0)
# [1] FALSE
res2 = drop.arc(res, "A", "B")
compare(res, res2, debug = TRUE)
# * nodes in the markov blanket of A in r1 not present in r2:
# [1] "B"
# * nodes in the markov blanket of A in r2 not present in r1:
# character(0)
# * nodes in the neighbourhood of A in r1 not present in r2:
# [1] "B"
# * nodes in the neighbourhood of A in r2 not present in r1:
# character(0)
# * nodes in the markov blanket of B in r1 not present in r2:
# [1] "A"
# * nodes in the markov blanket of B in r2 not present in r1:
# character(0)
# * nodes in the neighbourhood of B in r1 not present in r2:
# [1] "A"
# * nodes in the neighbourhood of B in r2 not present in r1:
# character(0)
# * undirected arcs in r1 not present in r2:
# [1] "A - B" "B - A"
# * undirected arcs in r2 not present in r1:
# character(0)
# [1] FALSE
res2 = reverse.arc(res, "A", "D")
compare(res, res2, debug = TRUE)
# * nodes in the markov blanket of A in r1 not present in r2:
# [1] "C"
# * nodes in the markov blanket of A in r2 not present in r1:
# character(0)
# * parents of A in r1 not present in r2:
# character(0)
# * parents of A in r2 not present in r1:
# [1] "D"
# * children of A in r1 not present in r2:
# [1] "D"
# * children of A in r2 not present in r1:
# character(0)
# * nodes in the markov blanket of C in r1 not present in r2:
# [1] "A"
# * nodes in the markov blanket of C in r2 not present in r1:
# character(0)
# * parents of D in r1 not present in r2:
# [1] "A"
# * parents of D in r2 not present in r1:
# character(0)
# * children of D in r1 not present in r2:
# character(0)
# * children of D in r2 not present in r1:
# [1] "A"
# * directed arcs in r1 not present in r2:
# [1] "A -> D"
# * directed arcs in r2 not present in r1:
# [1] "D -> A"
# [1] FALSE

## use debug = TRUE for more information
\dontrun{
set.arc(res, "A", "B", debug = TRUE)
drop.arc(res, "A", "B", debug = TRUE)
reverse.arc(res, "A", "D", debug = TRUE)
}
}
\author{ Marco Scutari }
\keyword{utilities}
\keyword{graphs}
