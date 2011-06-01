\name{arc.strength}
\alias{arc.strength}
\alias{boot.strength}
\alias{custom.strength}
\title{ Measure arc strength}
\description{

  Strength of the probabilistic relations expressed by the
  arcs of a Bayesian network.

}
\usage{
# strength of the arcs present in x.
arc.strength(x, data, criterion = NULL, ..., debug = FALSE)
# strength of all possible arcs, as learned from bootstrapped data.
boot.strength(data, R = 200, m = nrow(data),
  algorithm, algorithm.args = list(), cpdag = FALSE, debug = FALSE)
# strength of all possible arcs, from a list of custom networks. 
custom.strength(networks, nodes, cpdag = FALSE, debug = FALSE)
}
\arguments{
  \item{x}{an object of class \code{bn}.}
  \item{networks}{a list, containing either object of class \code{bn}
    or arc sets (matrices or data frames with two columns, optionally
    labeled "from" and "to").}
  \item{nodes}{a vector of character strings, the labels of the nodes
    in the network.}
  \item{data}{a data frame containing the data the Bayesian network
    was learned from.}
  \item{criterion}{a character string, the label of a score
    function, the label of an independence test or \code{bootstrap}.
    See \code{\link{bnlearn-package}} for details on the first two
    possibilities.}
  \item{R}{a positive integer, the number of bootstrap replicates.}
  \item{m}{a positive integer, the size of each bootstrap replicate.}
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

}
\value{

  \code{arc.strength}, \code{boot.strength} and \code{custom.strength}
  return an object of class \code{bn.strength}; \code{boot.strength}
  and \code{custom.strength} also include information about the relative
  probabilities of arc directions.
  See \code{\link{bn.strength class}} for details.

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
arcs[(arcs$strength > 0.85) & (arcs$direction > 0.5), ]
#    from to strength direction
# 1     A  B        1         1
# 3     A  D        1         1
# 9     B  E        1         1
# 13    C  D        1         1
# 30    F  E        1         1
\dontrun{
start = random.graph(nodes = names(learning.test), num = 50)
netlist = lapply(start, function(net) {
  hc(learning.test, score = "bde", iss = 10, start = net) })
arcs = custom.strength(netlist, nodes = names(learning.test))
arcs[(arcs$strength > 0.85) & (arcs$direction > 0.5), ]
#    from to strength direction
# 1     A  B        1      1.00
# 3     A  D        1      1.00
# 9     B  E        1      0.98
# 13    C  D        1      0.96
# 30    F  E        1      0.66
}
}
\seealso{\code{\link{strength.plot}}, \code{\link{choose.direction}},
  \code{\link{score}}, \code{\link{ci.test}}. }
\author{ Marco Scutari }
\keyword{htest}
