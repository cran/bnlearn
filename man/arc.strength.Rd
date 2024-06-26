\name{arc.strength}
\alias{arc.strength}
\alias{boot.strength}
\alias{custom.strength}
\alias{bf.strength}
\alias{mean.bn.strength}
\alias{averaged.network}
\alias{inclusion.threshold}
\title{Measure arc strength}
\description{

  Measure the strength of the probabilistic relationships expressed by the arcs
  of a Bayesian network, and use model averaging to build a network containing
  only the significant arcs.

}
\usage{
# strength of the arcs present in x.
arc.strength(x, data, criterion = NULL, ..., debug = FALSE)
# strength of all possible arcs, as learned from bootstrapped data.
boot.strength(data, cluster, R = 200, m = nrow(data),
  algorithm, algorithm.args = list(), cpdag = TRUE, shuffle = TRUE,
  debug = FALSE)
# strength of all possible arcs, from a list of custom networks.
custom.strength(networks, nodes, weights = NULL, cpdag = TRUE, debug = FALSE)
# strength of all possible arcs, computed using Bayes factors.
bf.strength(x, data, score, ..., debug = FALSE)

# average arc strengths.
\method{mean}{bn.strength}(x, ..., weights = NULL)

# averaged network structure.
averaged.network(strength, threshold)
# strength threshold for inclusion in the averaged network structure.
inclusion.threshold(strength)
}
\arguments{
  \item{x}{an object of class \code{bn.strength} (for \code{mean()}) or of
    class \code{bn} (for all other functions).}
  \item{networks}{a list, containing either object of class \code{bn} or arc
    sets (matrices or data frames with two columns, optionally labeled "from"
    and "to"); or an object of class \code{bn.kcv} or \code{bn.kcv.list}
    from \code{bn.cv()}.}
  \item{data}{a data frame containing the data the Bayesian network was
    learned from (for \code{arc.strength()}) or that will be used to compute
    the arc strengths (for \code{boot.strength()} and \code{bf.strength()}).}
  \item{cluster}{an optional cluster object from package \pkg{parallel}.}
  \item{strength}{an object of class \code{bn.strength}, see below.}
  \item{threshold}{a numeric value, the minimum strength required for an
    arc to be included in the averaged network. The default value is the
    \code{threshold} attribute of the \code{strength} argument.}
  \item{nodes}{a vector of character strings, the labels of the nodes in the
    network.}
  \item{criterion,score}{a character string. For \code{arc.strength()}, the
    label of a score function or an independence test; see
    \code{\link{network scores}} for details.} For \code{bf.strength()}, the
    label of the score used to compute the Bayes factors; see \code{\link{BF}}
    for details.
  \item{R}{a positive integer, the number of bootstrap replicates.}
  \item{m}{a positive integer, the size of each bootstrap replicate.}
  \item{weights}{a vector of non-negative numbers, to be used as weights
    when averaging arc strengths (in \code{mean()}) or network structures (in
    \code{custom.strength()}) to compute strength coefficients. If \code{NULL},
    weights are assumed to be uniform.}
  \item{cpdag}{a boolean value. If \code{TRUE} the (PDAG of) the equivalence
    class is used instead of the network structure itself. It should make it
    easier to identify score-equivalent arcs.}
  \item{shuffle}{a boolean value. If \code{TRUE} the columns in the data are
    permuted in each bootstrap sample to enforce the fact that the ordering of
    the variables in the data should be an invariant.}
  \item{algorithm}{a character string, the structure learning algorithm to be
    applied to the bootstrap replicates. See \code{\link{structure learning}}
    and the documentation of each algorithm for details.}
  \item{algorithm.args}{a list of extra arguments to be passed to the learning
    algorithm.}
  \item{\dots}{in \code{arc.strength()}, the additional tuning parameters for
    the network score (if \code{criterion} is the label of a score function,
    see \code{\link{score}} for details), the conditional independence test
    (currently the only one is \code{B}, the number of permutations). In
    \code{mean()}, additional objects of class \code{bn.strength} to average.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
}
\details{

  \code{arc.strength()} computes a measure of confidence or strength for each
  arc, while keeping fixed the rest of the network structure.

  If \code{criterion} is a conditional independence test, the strength is a
  p-value (so the lower the value, the stronger the relationship). The
  conditional independence test would be that to drop the arc from the network.
  The only two possible additional arguments are \code{alpha}, which sets the
  significance threshold that is used in \code{strength.plot()}; and \code{B},
  the number of permutations to be generated for each permutation test.

  If \code{criterion} is the label of a score function, the strength is
  measured by the score gain/loss which would be caused by the arc's removal.
  In other words, it is the difference between the score of the network
  in which the arc is not present and the score of the network in which the arc
  is present. Negative values correspond to decreases in the network score
  and positive values correspond to increases in the network score (the stronger
  the relationship, the more negative the difference). There may be additional
  arguments depending on the choice of the score, see \code{\link{score}} for
  details. The significance threshold is set to \code{0}.

  \code{boot.strength()} estimates the strength of each arc as its empirical
  frequency over a set of networks learned from bootstrap samples. It computes
  the probability of each arc (modulo its direction) and the probabilities of
  each arc's directions conditional on the arc being present in the graph (in
  either direction). The significance threshold is computed automatically from
  the strength estimates.

  \code{bf.strength()} estimates the strength of each arc using Bayes factors
  to overcome the fact that Bayesian posterior scores are not normalised, and
  uses the latter to estimate the probabilities of all possible states of an
  arc given the rest of the network. The significance threshold is set to
  \code{1}.

  \code{custom.strength()} takes a list of networks and estimates arc strength
  in the same way as \cr \code{boot.strength()}.

  Model averaging is supported for objects of class \code{bn.strength} returned
  by \code{\link{boot.strength}}, \code{\link{custom.strength}} and
  \code{\link{bf.strength}}. The returned network contains the arcs whose
  strength is greater than the \code{threshold} attribute of the
  \code{bn.strength} object passed to \code{averaged.network()}.

}
\note{

  \code{averaged.network()} typically returns a completely directed graph; an
  arc can be undirected if and only if the probability of each of its directions
  is exactly 0.5. This may happen, for example, if the arc is undirected in all
  the networks being averaged.

}
\value{

  \code{arc.strength()}, \code{boot.strength()}, \code{custom.strength()},
  \code{bf.strength()} and \code{mean()} return an object of class
  \code{bn.strength}; \code{boot.strength()} and \code{custom.strength()} also
  include information about the relative probabilities of arc directions.

  \code{averaged.network()} returns an object of class \code{bn}.

  See \code{\link{bn.strength class}} and \code{\link{bn-class}} for details.

}
\references{

  \bold{for model averaging and boostrap strength (confidence):}

  Friedman N, Goldszmidt M, Wyner A (1999). "Data Analysis with Bayesian
	  Networks: A Bootstrap Approach". \emph{Proceedings of the 15th Annual
	  Conference on Uncertainty in Artificial Intelligence}, 196--201.

  \bold{for the computation of the bootstrap strength (confidence) significance
    threshold:}

  Scutari M, Nagarajan R (2013). "On Identifying Significant Edges in Graphical
    Models of Molecular Networks". \emph{Artificial Intelligence in Medicine},
    \strong{57}(3):207--217.

}
\examples{
data(learning.test)
dag = hc(learning.test)
arc.strength(dag, learning.test)

\dontrun{
arcs = boot.strength(learning.test, algorithm = "hc")
arcs[(arcs$strength > 0.85) & (arcs$direction >= 0.5), ]
averaged.network(arcs)

start = random.graph(nodes = names(learning.test), num = 50)
netlist = lapply(start, function(net) {
  hc(learning.test, score = "bde", iss = 10, start = net) })
arcs = custom.strength(netlist, nodes = names(learning.test),
         cpdag = FALSE)
arcs[(arcs$strength > 0.85) & (arcs$direction >= 0.5), ]
modelstring(averaged.network(arcs))

bf.strength(dag, learning.test, score = "bds", prior = "marginal")
}}
\seealso{\code{\link{strength.plot}}, \code{\link{score}},
  \code{\link{ci.test}}.}
\author{Marco Scutari}
\keyword{network scores}
\keyword{independence tests}
\keyword{structure learning}
