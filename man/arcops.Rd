\name{arc operations}
\alias{arc operations}
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
  \item{x}{an object of class \code{bn}.}
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
      the arc \code{from} \eqn{\rightarrow}{->} \code{to} is added.
    \item if there is an undirected arc between \code{from} and
      \code{to}, its direction is set to \code{from}
      \eqn{\rightarrow}{-} \code{to}.
    \item if the arc \code{to} \eqn{\rightarrow}{->} \code{from}
      is present, it's reversed.
    \item if the arc \code{from} \eqn{\rightarrow}{->} \code{to}
      is present, no action is taken.
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
    \item if the arc \code{to} \eqn{\rightarrow}{->} \code{from}
      is present, it's reversed.
    \item if the arc \code{from} \eqn{\rightarrow}{->} \code{to}
      is present, it's reversed.
  }

}
\value{

  \code{set.arc} and \code{drop.arc} return invisibly an updated
  copy of \code{x}.

}
\examples{
data(learning.test)
res = gs(learning.test)

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
