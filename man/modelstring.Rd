\name{model string utilities}
\alias{model string utilities}
\alias{modelstring}
\alias{modelstring<-}
\alias{model2network}
\alias{as.character.bn}
\alias{as.bn}
\alias{as.bn.character}
\title{Build a model string from a Bayesian network and vice versa}
\description{

  Build a model string from a Bayesian network and vice versa.

}
\usage{
modelstring(x)
modelstring(x, debug = FALSE) <- value

model2network(string, ordering = NULL, debug = FALSE)

\method{as.character}{bn}(x, ...)
\method{as.bn}{character}(x)
}
\arguments{
  \item{x}{an object of class \code{bn}. \code{modelstring} (but not its
    replacement form) accepts also objects of class \code{bn.fit}.}
  \item{string}{a character string describing the Bayesian network.}
  \item{ordering}{the labels of all the nodes in the graph; their order is the
    node ordering used in the construction of the \code{bn} object. If
    \code{NULL} the nodes are sorted alphabetically.}
  \item{value}{a character string, the same as the \code{string}.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
  \item{\dots}{extra arguments from the generic method (currently ignored).}
}
\details{

  The strings returned by \code{modelstring} have the same format as the ones
  returned by the \code{modelstring} function in package \pkg{deal}; network
  structures may be easily exported to and imported from that package (via the
  \code{model2network} function).

  The format of the model strings is as follows. The local structure of each
  node is enclosed in square brackets ("\code{[]}"); the first string is the
  label of that node. The parents of the node (if any) are listed after a
  ("\code{|}") and separated by colons ("\code{:}"). All nodes (including
  isolated and root nodes) must be listed.

}
\value{

  \code{model2network} and \code{as.bn} return an object of class \code{bn};
  \code{modelstring} and \code{as.character.bn} return a character string.

}
\examples{
data(learning.test)
res = set.arc(gs(learning.test), "A", "B")
res
modelstring(res)
res2 = model2network(modelstring(res))
res2
all.equal(res, res2)
}
\author{ Marco Scutari }
\keyword{utilities}
\keyword{graphs}
