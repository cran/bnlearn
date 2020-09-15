\name{constraint-based algorithms}
\alias{constraint-based algorithms}
\alias{pc.stable}
\alias{gs}
\alias{iamb}
\alias{fast.iamb}
\alias{inter.iamb}
\alias{iamb.fdr}
\alias{mmpc}
\alias{si.hiton.pc}
\alias{hpc}
\title{Constraint-based structure learning algorithms}
\description{

  Learn the equivalence class of a directed acyclic graph (DAG) from data using
  the PC, Grow-Shrink (GS), Incremental Association (IAMB), Fast Incremental
  Association (Fast-IAMB), Interleaved Incremental Association (Inter-IAMB),
  Incremental Association with FDR (IAMB-FDR), Max-Min Parents and Children
  (MMPC), Semi-Interleaved HITON-PC or Hybrid Parents and Children (HPC)
  constraint-based algorithms.
}
\usage{
pc.stable(x, cluster = NULL, whitelist = NULL, blacklist = NULL, test = NULL,
  alpha = 0.05, B = NULL, max.sx = NULL, debug = FALSE, undirected = FALSE)
gs(x, cluster = NULL, whitelist = NULL, blacklist = NULL, test = NULL,
  alpha = 0.05, B = NULL, max.sx = NULL, debug = FALSE, undirected = FALSE)
iamb(x, cluster = NULL, whitelist = NULL, blacklist = NULL, test = NULL,
  alpha = 0.05, B = NULL, max.sx = NULL, debug = FALSE, undirected = FALSE)
fast.iamb(x, cluster = NULL, whitelist = NULL, blacklist = NULL, test = NULL,
  alpha = 0.05, B = NULL, max.sx = NULL, debug = FALSE, undirected = FALSE)
inter.iamb(x, cluster = NULL, whitelist = NULL, blacklist = NULL, test = NULL,
  alpha = 0.05, B = NULL, max.sx = NULL, debug = FALSE, undirected = FALSE)
iamb.fdr(x, cluster = NULL, whitelist = NULL, blacklist = NULL, test = NULL,
  alpha = 0.05, B = NULL, max.sx = NULL, debug = FALSE, undirected = FALSE)
mmpc(x, cluster = NULL, whitelist = NULL, blacklist = NULL, test = NULL,
  alpha = 0.05, B = NULL, max.sx = NULL, debug = FALSE, undirected = TRUE)
si.hiton.pc(x, cluster = NULL, whitelist = NULL, blacklist = NULL, test = NULL,
  alpha = 0.05, B = NULL, max.sx = NULL, debug = FALSE, undirected = TRUE)
hpc(x, cluster = NULL, whitelist = NULL, blacklist = NULL, test = NULL,
  alpha = 0.05, B = NULL, max.sx = NULL, debug = FALSE, undirected = TRUE)
}
\arguments{
  \item{x}{a data frame containing the variables in the model.}
  \item{cluster}{an optional cluster object from package \pkg{parallel}.}
  \item{whitelist}{a data frame with two columns (optionally labeled "from"
    and "to"), containing a set of arcs to be included in the graph.}
  \item{blacklist}{a data frame with two columns (optionally labeled "from"
    and "to"), containing a set of arcs not to be included in the graph.}
  \item{test}{a character string, the label of the conditional independence
    test to be used in the algorithm. If none is specified, the default test
    statistic is the \emph{mutual information} for categorical variables, the
    Jonckheere-Terpstra test for ordered factors and the \emph{linear
    correlation} for continuous variables. See \code{\link{independence tests}}
    for details.}
  \item{alpha}{a numeric value, the target nominal type I error rate.}
  \item{B}{a positive integer, the number of permutations considered for each
    permutation test. It will be ignored with a warning if the conditional
    independence test specified by the \code{test} argument is not a
    permutation test.}
  \item{max.sx}{a positive integer, the maximum allowed size of the conditioning
    sets used in conditional independence tests. The default is that there is
    no limit on size.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
  \item{undirected}{a boolean value. If \code{TRUE} no attempt will be made to
    determine the orientation of the arcs; the returned (undirected) graph
    will represent the underlying structure of the Bayesian network.}
}
\note{

  Note that even when \code{undirected} is set to \code{FALSE} there is no
  guarantee that all arcs in the returned network will be directed; some arc
  directions are impossible to learn just from data due to score equivalence.
  \code{cextend()} provides a consistent extension of partially directed
  networks into directed acyclic graphs, which can then be used (for instance)
  for parameter learning.

  See \code{\link{structure learning}} for a complete list of structure learning
  algorithms with the respective references.

}
\value{

  An object of class \code{bn}.
  See \code{\link{bn-class}} for details.

}
\author{Marco Scutari}
\seealso{\link{local discovery algorithms},
  \link{score-based algorithms}, \link{hybrid algorithms}, \link{cextend}.}
\keyword{structure learning}
