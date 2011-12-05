\name{arc.strength}
\alias{arc.strength}
\alias{boot.strength}
\alias{custom.strength}
\alias{averaged.network}
\title{Measure arc strength}
\description{

  Measure the strength of the probabilistic relationships expressed
  by the arcs of a Bayesian network, and use model averaging to build
  a network containing only the significant arcs.

}
\usage{
# strength of the arcs present in x.
arc.strength(x, data, criterion = NULL, ..., debug = FALSE)
# strength of all possible arcs, as learned from bootstrapped data.
boot.strength(data, R = 200, m = nrow(data),
  algorithm, algorithm.args = list(), cpdag = TRUE, debug = FALSE)
# strength of all possible arcs, from a list of custom networks.
custom.strength(networks, nodes, weights = NULL, cpdag = TRUE,
  debug = FALSE)

# averaged network structure.
averaged.network(strength, nodes, threshold)
}
\arguments{
  \item{x}{an object of class \code{bn}.}
  \item{networks}{a list, containing either object of class \code{bn}
    or arc sets (matrices or data frames with two columns, optionally
    labeled "from" and "to").}
  \item{data}{a data frame containing the data the Bayesian network
    was learned from.}
  \item{strength}{an object of class \code{bn.strength}, see below.}
  \item{threshold}{a numeric value, the minimum strength required for
    an arc to be included in the averaged network. The default value
    is the \code{threshold} attribute of the \code{strength} argument.}
  \item{nodes}{a vector of character strings, the labels of the nodes
    in the network. In \code{\link{averaged.network}}, it defaults to
    the set of the unique node labels in the \code{strength} argument.}
  \item{criterion}{a character string, the label of a score
    function, the label of an independence test or \code{bootstrap}.
    See \code{\link{bnlearn-package}} for details on the first two
    possibilities.}
  \item{R}{a positive integer, the number of bootstrap replicates.}
  \item{m}{a positive integer, the size of each bootstrap replicate.}
  \item{weights}{a vector of non-negative numbers, to be used as
    weights when averaging network structures to compute strength
    coefficients. If \code{NULL}, weights are assumed to be
    uniform.}
  \item{cpdag}{a boolean value. If \code{TRUE} the (PDAG of) the
    equivalence class is used instead of the network structure
    itself. It should make it easier to identify score-equivalent
    arcs.}
  \item{algorithm}{a character string, the learning algorithm to be
    applied to the bootstrap replicates. Possible values are \code{gs},
    \code{iamb}, \code{fast.iamb}, \code{inter.iamb}, \code{mmpc},
    \code{hc}, \code{tabu}, \code{mmhc} and \code{rsmax2}. See
    \code{\link{bnlearn-package}} and the documentation of each algorithm
    for details.}
  \item{algorithm.args}{a list of extra arguments to be passed to
    the learning algorithm.}
  \item{\dots}{additional tuning parameters for the network score (if
    \code{criterion} is the label of a score function, see
    \code{\link{score}} for details), the conditional independence test
    (currently the only one is \code{B}, the number of permutations)
    or the bootstrap simulation (if \code{criterion} is set to
    \code{bootstrap}, see \code{\link{boot.strength}} for details).}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output
    is printed; otherwise the function is completely silent.}
}
\details{

  If \code{criterion} is a conditional independence test, the strength
  is a p-value (so the lower the value, the stronger the relationship).
  The only possible additional parameter is \code{B}, the number of
  permutations to be generated for each permutation test.

  If \code{criterion} is the label of a score function, the strength is
  measured by the score gain/loss which would be caused by the arc's
  removal. There may be additional parameters depending on the choice
  of the score, see \code{\link{score}} for details.

  If \code{criterion} is \code{bootstrap}, the strength is computed as
  in \code{\link{boot.strength}}. The additional parameters are \code{R},
  \code{m}, \code{algorithm} and \code{algorithm.args}; if the latter two
  are not specified, the values stored in \code{x} are used.

  Model averaging is supported for objects of class \code{bn.strength}
  returned by \code{\link{boot.strength}}, by \code{\link{custom.strength}},
  or by \code{\link{arc.strength}} with \code{criterion} set to
  \code{bootstrap}. The returned network contains the arcs whose strength
  is greater than the \code{threshold} attribute of the \code{bn.strength}
  object passed to \code{averaged.network}.

}
\value{

  \code{arc.strength}, \code{boot.strength} and \code{custom.strength}
  return an object of class \code{bn.strength}; \code{boot.strength}
  and \code{custom.strength} also include information about the relative
  probabilities of arc directions.

  \code{averaged.network} returns an object of class \code{bn}.

  See \code{\link{bn.strength class}} and \code{\link{bn-class}} for details.

}
\references{

  \bold{for model averaging and boostrap strength (confidence):}

  Friedman N, Goldszmidt M, Wyner A (1999). "Data Analysis with Bayesian
      Networks: A Bootstrap Approach". In "UAI '99: Proceedings of the
      15th Annual Conference on Uncertainty in Artificial Intelligence",
      pp. 196-20. Morgan Kaufmann.

  \bold{for the computation of the strength (confidence) significance threshold:}

  Scutari M, Nagarajan R (2011). "On Identifying Significant Edges in
      Graphical Models". In "Proceedings of the Workshop 'Probabilistic
      Problem Solving in Biomedicine' of the 13th Artificial Intelligence
      in Medicine (AIME) Conference", pp. 15-27.

}
\examples{
data(learning.test)
res = gs(learning.test)
res = set.arc(res, "A", "B")
arc.strength(res, learning.test)
#   from to      strength
# 1    A  B  0.000000e+00
# 2    A  D  0.000000e+00
# 3    B  E 1.024198e-320
# 4    C  D  0.000000e+00
# 5    F  E 3.935648e-245
arcs = boot.strength(learning.test, algorithm = "hc")
arcs[(arcs$strength > 0.85) & (arcs$direction >= 0.5), ]
#    from to strength direction
# 1     A  B        1       0.5
# 3     A  D        1       1.0
# 6     B  A        1       0.5
# 9     B  E        1       1.0
# 13    C  D        1       1.0
# 30    F  E        1       1.0
averaged.network(arcs)
#
#   Random/Generated Bayesian network
#
#   model:
#    [A][C][F][B|A][D|A:C][E|B:F]
#   nodes:                                 6
#   arcs:                                  5
#     undirected arcs:                     0
#     directed arcs:                       5
#   average markov blanket size:           2.33
#   average neighbourhood size:            1.67
#   average branching factor:              0.83
#
#   generation algorithm:                  Model Averaging
#   significance threshold:                0.025

\dontrun{
start = random.graph(nodes = names(learning.test), num = 50)
netlist = lapply(start, function(net) {
  hc(learning.test, score = "bde", iss = 10, start = net) })
arcs = custom.strength(netlist, nodes = names(learning.test))
arcs[(arcs$strength > 0.85) & (arcs$direction >= 0.5), ]
#    from to strength direction
# 1     A  B        1      1.00
# 3     A  D        1      1.00
# 9     B  E        1      0.98
# 13    C  D        1      0.96
# 30    F  E        1      0.66
modelstring(averaged.network(arcs))
# [1] "[A][C][F][B|A][D|A:C][E|B:F]"
}
}
\seealso{\code{\link{strength.plot}}, \code{\link{choose.direction}},
  \code{\link{score}}, \code{\link{ci.test}}. }
\author{ Marco Scutari }
\keyword{htest}
