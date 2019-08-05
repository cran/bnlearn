\name{score-based algorithms}
\alias{score-based algorithms}
\alias{hc}
\alias{tabu}
\title{Score-based structure learning algorithms}
\description{

  Learn the structure of a Bayesian network using a hill-climbing (HC) or a
  Tabu search (TABU) greedy search.

}
\usage{
hc(x, start = NULL, whitelist = NULL, blacklist = NULL, score = NULL, ...,
  debug = FALSE, restart = 0, perturb = 1, max.iter = Inf, maxp = Inf, optimized = TRUE)
tabu(x, start = NULL, whitelist = NULL, blacklist = NULL, score = NULL, ...,
  debug = FALSE, tabu = 10, max.tabu = tabu, max.iter = Inf, maxp = Inf, optimized = TRUE)
}
\arguments{
  \item{x}{a data frame containing the variables in the model.}
  \item{start}{an object of class \code{bn}, the preseeded directed acyclic
    graph used to initialize the algorithm. If none is specified, an empty one
    (i.e. without any arc) is used.}
  \item{whitelist}{a data frame with two columns (optionally labeled "from" and
    "to"), containing a set of arcs to be included in the graph.}
  \item{blacklist}{a data frame with two columns (optionally labeled "from" and
    "to"), containing a set of arcs not to be included in the graph.}
  \item{score}{a character string, the label of the network score to be used in
    the algorithm. If none is specified, the default score is the \emph{Bayesian
    Information Criterion} for both discrete and continuous data sets. See
    \code{\link{network scores}} for details.}
  \item{\dots}{additional tuning parameters for the network score. See
    \code{\link{score}} for details.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
  \item{restart}{an integer, the number of random restarts.}
  \item{tabu}{a positive integer number, the length of the tabu list used in the
    \code{tabu} function.}
  \item{max.tabu}{a positive integer number, the iterations tabu search can
    perform without improving the best network score.}
  \item{perturb}{an integer, the number of attempts to randomly
    insert/remove/reverse an arc on every random restart.}
  \item{max.iter}{an integer, the maximum number of iterations.}
  \item{maxp}{the maximum number of parents for a node. The default value is
    \code{Inf}.}
  \item{optimized}{a boolean value. If \code{TRUE} (the default), score caching
    is used to speed up structure learning.}
}
\note{

  See \code{\link{structure learning}} for a complete list of structure learning
  algorithms with the respective references.

}
\value{

  An object of class \code{bn}. See \code{\link{bn-class}} for details.

}
\author{Marco Scutari}
\seealso{\link{constraint-based algorithms}, \link{hybrid algorithms}, \cr
  \link{local discovery algorithms}, \link{alpha.star}.}
\keyword{structure learning}
