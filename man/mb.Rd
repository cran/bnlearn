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
\alias{amat}
\alias{amat<-}
\alias{alst}
\alias{alst<-}
\alias{parents}
\alias{parents<-}
\alias{children}
\alias{children<-}
\alias{spouses}
\alias{ancestors}
\alias{descendants}
\alias{in.degree}
\alias{out.degree}
\alias{degree}
\alias{degree,bn-method}
\alias{degree,bn.fit-method}
\alias{degree,bn.naive-method}
\alias{degree,bn.tan-method}
\alias{root.nodes}
\alias{leaf.nodes}
\alias{isolated.nodes}
\alias{nparams}
\alias{ntests}
\title{Miscellaneous utilities}
\description{

  Assign or extract various quantities of interest from an object of class
  \code{bn} of \code{bn.fit}.

}
\usage{
## nodes
mb(x, node)
nbr(x, node, max.dist = 1)
parents(x, node)
parents(x, node, debug = FALSE) <- value
children(x, node)
children(x, node, debug = FALSE) <- value
spouses(x, node)
ancestors(x, node)
descendants(x, node)
in.degree(x, node)
out.degree(x, node)
root.nodes(x)
leaf.nodes(x)
isolated.nodes(x)
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

## arc list
alst(x)
alst(x, check.cycles = TRUE, check.illegal = TRUE, debug = FALSE) <- value

## graphs
nparams(x, data, estimator = NULL, debug = FALSE)
ntests(x)

## shared with the graph package.
# these used to be a simple nodes(x) function.
\S4method{nodes}{bn}(object)
\S4method{nodes}{bn.fit}(object)
# these used to be a simple degree(x, node) function.
\S4method{degree}{bn}(object, Nodes)
\S4method{degree}{bn.fit}(object, Nodes)
}
\arguments{
  \item{x,object}{an object of class \code{bn} or \code{bn.fit}. The replacement
    form of \code{parents()}, \code{children()}, \code{arcs()}, \code{amat} and
    \code{alst()} requires an object of class \code{bn}.}
  \item{node,Nodes}{a character string, the label of a node.}
  \item{max.dist}{a positive integer, the maximum number of arcs (distance) from
    the node at the center of the neighbourhood.}
  \item{value}{either a vector of character strings (for \code{parents()} and
    \code{children()}), an adjacency matrix (for \code{amat()}), an arc list
    (for \code{alst()}) or a data frame with two columns (optionally labelled
    "from" and "to", for \code{arcs()}).}
  \item{data}{a data frame containing the data the Bayesian network was learned
    from. It's only needed if \code{x} is an object of class \code{bn}.}
  \item{estimator}{a character string, the parameter estimator that will be used
    to fit the network; see \code{\link{bn.fit}} for the available estimators. It
    determines the number of parameters of the continuous nodes. If \code{NULL},
    it defaults to the estimator implied by the type of the data.}
  \item{check.cycles}{a boolean value. If \code{FALSE}, the returned network
    will not be checked for cycles.}
  \item{check.illegal}{a boolean value. If \code{TRUE}, arcs that break the
    parametric assumptions of \code{x}, such as those from continuous to
    discrete nodes in conditional Gaussian networks, cause an error.}
  \item{debug}{a boolean value. If \code{TRUE}, a lot of debugging output is
    printed. Otherwise, the function is completely silent.}
}
\details{

  \itemize{

    \item The number of parameters of a discrete Bayesian network is defined as
      the sum of the number of logically independent parameters of each node
      given its parents (Chickering, 1995).
    \item For Gaussian Bayesian networks, the distribution of each node can be
      viewed as a linear regression, so it has a number of parameters equal to
      the number of the node's parents plus one (the intercept), as per
      Neapolitan (2003).
   \item For conditional linear Gaussian networks, the number of parameters of
     discrete and Gaussian nodes is as above. The number of parameters of
     conditional Gaussian nodes is equal to \code{1} plus the number of
     continuous parents (who get one regression coefficient each, plus the
     intercept) times the number of configurations of the discrete parents (each
     configuration has an associated regression model).
   \item For zero-inflated networks, the number of parameters of each node is
     equal to twice the number of parents plus three. They correspond to the
     set of regression coefficients of the zero-inflation component, that of
     the counts component, and the shape parameter (the dispersion in the
     hyper-Poisson, the number of failures in the negative binomial).

  }

  Both Gaussian and zero-inflated networks are learned from numeric data, so the
  number of parameters of the continuous nodes is ambiguous. The
  \code{estimator} argument selects which parametric assumption is used to count
  them; by default it is the estimator implied by the type of the data (which
  amounts to the Gaussian parametrisation for numeric data).

  \pkg{bnlearn} implements arc lists (equivalently adjacency lists) as a list
  with one element per node. Elements are named after the nodes in \code{x}, and
  contain a character vector holding the labels of the children of that node.
  Undirected arcs are treated as if they were bi-directed, so each node is
  listed as the other's child.

}
\value{

  \code{mb()}, \code{nbr()}, \code{nodes()}, \code{parents()},
  \code{children()}, \code{spouses()}, \code{ancestors()}, \code{descendants()},
  \code{root.nodes()}, \code{leaf.nodes()} and \code{isolated.nodes()} return a
  vector of character strings.

  \code{arcs()}, \code{directed.arcs()}, \code{undirected.arcs()},
  \code{incoming.arcs()}, \code{outgoing.arcs()}, \code{incident.arcs()}, \cr
  \code{compelled.arcs()}, \code{reversible.arcs()}, return a matrix of two
  columns of character strings.

  \code{narcs()} and \code{nnodes()} return the number of arcs and nodes in the
  graph, respectively.

  \code{amat()} returns a matrix of 0/1 integer values.

  \code{alst()} returns a list of character vectors.

  \code{degree()}, \code{in.degree()}, \code{out.degree()}, \code{nparams()} and
  \code{ntests()} return an integer.

}
\references{

  Chickering DM (1995). "A Transformational Characterization of Equivalent
    Bayesian Network Structures." \emph{Proceedings of the Eleventh Annual
    Conference on Uncertainty in Artificial Intelligence}, 87--98.

  Neapolitan RE (2003). \emph{Learning Bayesian Networks}. Prentice Hall.

}
\examples{
data(learning.test)
cpdag = pc.stable(learning.test)

##  the Markov blanket of A.
mb(cpdag, "A")
## the neighbourhood of F.
nbr(cpdag, "F")
## the arcs in the graph.
arcs(cpdag)
## the nodes of the graph.
nodes(cpdag)
## the adjacency matrix for the nodes of the graph.
amat(cpdag)
## the adjacency list of a graph.
alst(cpdag)
## the parents of D.
parents(cpdag, "D")
## the children of A.
children(cpdag, "A")
## the root nodes of the graph.
root.nodes(cpdag)
## the leaf nodes of the graph.
leaf.nodes(cpdag)
## number of parameters of the Bayesian network.
dag = set.arc(cpdag, "A", "B")
nparams(dag, learning.test)
}
\author{Marco Scutari}
\keyword{convenience functions}
\keyword{graphs}
