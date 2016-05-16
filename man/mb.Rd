\name{misc utilities}
\alias{misc utilities}
\alias{mb}
\alias{nbr}
\alias{arcs}
\alias{arcs<-}
\alias{directed.arcs}
\alias{undirected.arcs}
\alias{incoming.arcs}
\alias{outgoing.arcs}
\alias{incident.arcs}
\alias{compelled.arcs}
\alias{reversible.arcs}
\alias{narcs}
\alias{nnodes}
\alias{nodes}
\alias{nodes,bn-method}
\alias{nodes,bn.fit-method}
\alias{nodes,bn.naive-method}
\alias{nodes,bn.tan-method}
\alias{nodes<-}
\alias{nodes<-,bn-method}
\alias{nodes<-,bn.fit-method}
\alias{nodes<-,bn.naive-method}
\alias{nodes<-,bn.tan-method}
\alias{amat}
\alias{amat<-}
\alias{parents}
\alias{parents<-}
\alias{children}
\alias{children<-}
\alias{in.degree}
\alias{out.degree}
\alias{degree}
\alias{degree,bn-method}
\alias{degree,bn.fit-method}
\alias{degree,bn.naive-method}
\alias{degree,bn.tan-method}
\alias{root.nodes}
\alias{leaf.nodes}
\alias{nparams}
\alias{ntests}
\alias{whitelist}
\alias{blacklist}
\title{Miscellaneous utilities}
\description{

  Assign or extract various quantities of interest from an object of class
  \code{bn} of \code{bn.fit}.

}
\usage{
## nodes
mb(x, node)
nbr(x, node)
parents(x, node)
parents(x, node, debug = FALSE) <- value
children(x, node)
children(x, node, debug = FALSE) <- value
in.degree(x, node)
out.degree(x, node)
root.nodes(x)
leaf.nodes(x)
nnodes(x)

## arcs
arcs(x)
arcs(x, check.cycles = TRUE, check.illegal = TRUE, debug = FALSE) <- value
directed.arcs(x)
undirected.arcs(x)
incoming.arcs(x, node)
outgoing.arcs(x, node)
incident.arcs(x, node)
compelled.arcs(x)
reversible.arcs(x)
narcs(x)

## adjacency matrix
amat(x)
amat(x, check.cycles = TRUE, check.illegal = TRUE, debug = FALSE) <- value

## graphs
nparams(x, data, effective = FALSE, debug = FALSE)
ntests(x)
whitelist(x)
blacklist(x)

## shared with the graph package.
# these used to be a simple nodes(x) function.
\S4method{nodes}{bn}(object)
\S4method{nodes}{bn.fit}(object)
# these used to be a simple degree(x, node) function.
\S4method{degree}{bn}(object, Nodes)
\S4method{degree}{bn.fit}(object, Nodes)
# re-label the nodes.
\S4method{nodes}{bn}(object) <- value
\S4method{nodes}{bn.fit}(object) <- value
}
\arguments{
  \item{x,object}{an object of class \code{bn} or \code{bn.fit}. The replacement
    form of \code{parents}, \code{children}, \code{arcs} and \code{amat}
    requires an object of class \code{bn}.}
  \item{node,Nodes}{a character string, the label of a node.}
  \item{value}{either a vector of character strings (for \code{parents} and
    \code{children}), an adjacency matrix (for \code{amat}) or a data frame with
    two columns (optionally labeled "from" and "to", for \code{arcs}).}
  \item{data}{a data frame containing the data the Bayesian network was learned
    from. It's only needed if \code{x} is an object of class \code{bn}.}
  \item{check.cycles}{a boolean value. If \code{FALSE} the returned network will
    not be checked for cycles.}
  \item{check.illegal}{a boolean value. If \code{TRUE} arcs that break the 
    parametric assumptions of \code{x}, such as those from continuous to
    discrete nodes in conditional Gaussian networks, cause an error.}
  \item{effective}{a boolean value. If \code{TRUE} the number of non-zero free
    parameters is returned, that is, the effective degrees of freedom of the
    network; otherwise the theoretical number of parameters is returned.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
}
\details{

  The number of parameters of a discrete Bayesian network is defined as the
  sum of the number of logically independent parameters of each node given its
  parents (Chickering, 1995). For Gaussian Bayesian networks the distribution
  of each node can be viewed as a linear regression, so it has a number of
  parameters equal to the number of the parents of the node plus one (the
  intercept) as per Neapolitan (2003).

}
\value{

  \code{mb}, \code{nbr}, \code{nodes}, \code{parents}, \code{children},
  \code{root.nodes} and \code{leaf.nodes} return a vector of character strings.

  \code{arcs}, \code{directed.arcs}, \code{undirected.arcs},
  \code{incoming.arcs}, \code{outgoing.arcs}, \code{incident.arcs}, \cr
  \code{compelled.arcs}, \code{reversible.arcs}, \code{whitelist} and
  \code{blacklist} return a matrix of two columns of character strings.

  \code{narcs} and \code{nnodes} return the number of arcs and nodes in the
  graph, respectively.

  \code{amat} returns a matrix of 0/1 integer values.

  \code{degree}, \code{in.degree}, \code{out.degree}, \code{nparams} and
  \code{ntests} return an integer.

}
\references{

  Chickering DM (1995). "A Transformational Characterization of Equivalent
      Bayesian Network Structures". In "UAI '95: Proceedings of the Eleventh
      Annual Conference on Uncertainty in Artificial Intelligence", pp. 87-98.
      Morgan Kaufmann.

  Neapolitan RE (2003). \emph{Learning Bayesian Networks}. Prentice Hall.

}
\examples{
data(learning.test)
res = gs(learning.test)

##  the Markov blanket of A.
mb(res, "A")
## the neighbourhood of F.
nbr(res, "F")
## the arcs in the graph.
arcs(res)
## the nodes of the graph.
nodes(res)
## the adjacency matrix for the nodes of the graph.
amat(res)
## the parents of D.
parents(res, "D")
## the children of A.
children(res, "A")
## the root nodes of the graph.
root.nodes(res)
## the leaf nodes of the graph.
leaf.nodes(res)
## number of parameters of the Bayesian network.
res = set.arc(res, "A", "B")
nparams(res, learning.test)
}
\author{ Marco Scutari }
\keyword{utilities}
\keyword{graphs}
