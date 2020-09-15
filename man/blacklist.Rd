\name{utilities for whitelists and blacklists}
\alias{whitelist}
\alias{blacklist}
\alias{ordering2blacklist}
\alias{tiers2blacklist}
\alias{set2blacklist}
\title{Get or create whitelists and blacklists}
\description{

  Extract whitelists and blacklists from an object of class \code{bn}, or create
  them for use in structure learning.

}
\usage{
whitelist(x)
blacklist(x)

ordering2blacklist(nodes)
tiers2blacklist(tiers)
set2blacklist(set)
}
\arguments{
  \item{x}{an object of class \code{bn}.}
  \item{nodes,set}{a vector of character strings, the labels of the nodes.}
  \item{tiers}{a vector of character strings or a list, see below.}

}
\details{

  \code{ordering2blacklist()} takes a vector of character strings (the labels
  of the nodes), which specifies a complete node ordering. An object of class
  \code{bn} or \code{bn.fit}; in that case, the node ordering is derived by the
  graph. In both cases, the blacklist returned by \code{ordering2blacklist()}
  contains all the possible arcs that violate the specified node ordering.

  \code{tiers2blacklist()} takes (again) a vector of character strings (the
  labels of the nodes), which specifies a complete node ordering, or a list of
  character vectors, which specifies a partial node ordering. In the latter
  case, all arcs going from a node in a particular element of the list
  (sometimes known as \emph{tier}) to a node in one of the previous elements
  are blacklisted. Arcs between nodes in the same element are not blacklisted.

  \code{set2blacklist()} creates a blacklist containing all the arcs between any
  two of the nodes whose labels are passed as the argument \code{set}.

}
\value{

  \code{whitelist()} and \code{blacklist()} return a matrix of character string
  with two columns, named \code{from} and \code{to}, if whitelist or a blacklist
  have been used to learn the \code{bn} object passed as their argument.

  \code{ordering2blacklist()}, \code{tiers2blacklist()} and
  \code{set2blacklist()} return a sanitized \code{blacklist} (a two-column
  matrix, whose columns are labeled \code{from} and \code{to}).

}
\examples{
tiers2blacklist(list(LETTERS[1:3], LETTERS[4:6]))
set2blacklist(LETTERS[1:3])
ordering2blacklist(LETTERS[1:6])
}

\author{Marco Scutari}
\keyword{convenience functions}
\keyword{structure learning}
