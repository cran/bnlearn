\name{bnlearn-package}
\alias{bnlearn-package}
\alias{bnlearn}
\docType{package}
\title{

  Bayesian network constraint-based structure learning.

}
\description{

  Bayesian network learning via constraint-based algorithms.

}
\details{

  \tabular{ll}{

    Package: \tab bnlearn\cr
    Type: \tab Package\cr
    Version: \tab 0.7\cr
    Date: \tab 2008-03-25\cr
    License: \tab GPLv2 or later

  }

  This package implements some \emph{constraint-based} algorithms for learning
  the structure of Bayesian networks. Also known as \emph{conditional
  independence learners}, they are all optimized derivatives of the
  \emph{Inductive Causation} algorithm (Verma and Pearl, 1991).

  These algorithms differ in the way they detect the Markov blankets of
  the variables, which in turn are used to compute the structure of the
  Bayesian network. Proofs of correctness are present in the respective
  papers.

  A \emph{score-based} learning algorithm (greedy search via hill-climbing) is
  implemented as well for comparison purposes.

}
\author{

  Marco Scutari\cr
  Department of Statistical Sciences\cr
  University of Padua

  Maintainer: Marco Scutari \email{marco.scutari@gmail.com}

}
\section{Available constraint-based learning algorithms}{

  \itemize{

    \item \emph{Grow-Shrink} (\code{\link{gs}}): based on the \emph{Grow-Shrink
      Markov Blanket} algorithm, a forward-selection technique for Markov blanket
      detection (Margaritis, 2003).
    \item \emph{Incremental Association} (\code{\link{iamb}}): based on the
      Markov blanket detection algorithm of the same name, which uses stepwise
      selection (Tsamardinos et al., 2003).
    \item \emph{Fast Incremental Association} (\code{\link{fast.iamb}}): a
      variant of IAMB which uses speculative forward-selection to reduce the
      number of conditional independence tests (Yaramakala and Margaritis, 2005).
    \item \emph{Interleaved Incremental Association} (\code{\link{inter.iamb}}):
      another variant of IAMB which interleaves the backward selection with
      the forward one to avid false positives (i.e. nodes erroneously included
      in the Markov blanket) (Yaramakala and Margaritis, 2005).

  }

  This package includes three implementations of each algorithm:

  \itemize{

    \item an optimized implementation (used when the \code{optimized}
      parameter is set to \code{TRUE}), which uses backtracking to roughly
      halve the number of independence tests.
    \item an unoptimized implementation (used when the \code{optimized}
      parameter is set to \code{FALSE}) which is better at uncovering
      possible erratic behaviour of the statistical tests.
    \item a cluster-aware implementation, which requires a running cluster
      set up with the \code{makeCluster} function from the \pkg{snow}
      package. See \link{snow integration} for a sample usage.

  }

  Their computational complexity of this class of algorithms is polynomial
  in the number of tests, usually \eqn{O(N^2)}{O(N^2)} (\eqn{O(N^4)}{O(N^4)}
  in the worst case scenario). The execution time also scales linearly with
  the size of the data set.

}
\section{Other learning algorithms}{

  \itemize{
    \item \emph{Hill-Climbing} (\code{\link{hc}}): an \emph{hill climbing}
      greedy search on the space of the directed graphs. The optimized
      implementation uses score caching, score decomposability and score
      equivalence to reduce the number of duplicated tests.}

}
\section{Available (conditional) independence tests}{

  The conditional independence tests used in \emph{constraint-based}
  algorithms in practice are statistical tests on the data set. Available
  tests (and the respective labels) are:

  \itemize{

    \item \emph{categorical data} (multinomial distribution)
      \itemize{

        \item \emph{mutual information} (\code{mi}): an information-theoretic
          distance measure. Is proportional to the log-likelihood ratio
          (they differ by a \eqn{2n}{2n} factor) and is related to the
          deviance of the tested models.
        \item \emph{fast mutual information} (\code{fmi}): a variant of the
          mutual information which is set to zero when there aren't at least
          five data per parameter.
        \item \emph{Akaike Information Criterion} (\code{aict}): an AIC-based
          indepence test, computed comparing the mutual information
          and the expected information gain.
        \item \emph{Cochran-Mantel-Haenszel} (\code{mh}): a stratified
          independence test, included for testing purposes only. See
          \code{\link{mantelhaen.test}} in package \code{\link{stats}}.

      }
    \item \emph{numerical data} (multivariate normal distribution)
      \itemize{

        \item \emph{linear correlation} (\code{cor}): linear correlation.
        \item \emph{Fisher's Z} (\code{zf}): a transformation of the
           linear correlation with asymptotic normal distribution. Used
           by commercial software (such as TETRAD II) for the PC algorithm
           (an R implementation is present in the \code{pcalg} package
           on CRAN).

      }

  }

}
\section{Available network scores}{

  Available scores (and the respective labels) are:

  \itemize{

    \item \emph{categorical data} (multinomial distribution)
      \itemize{

        \item the multinomial \emph{likelihood} (\code{lik}) score.
        \item the multinomial \emph{loglikelihood} (\code{loglik}) score.
        \item the \emph{Akaike Information Criterion} score (\code{aic}).
        \item the \emph{Bayesian Information Criterion} score (\code{bic}).
        \item the logarithm of the \emph{Bayesian Dirichlet equivalent}
          score (\code{bde} o \code{dir}), a score-equivalent dirichlet
          posterior density.
        \item the logarithm of the \emph{K2} score (\code{k2}), a dirichlet
          posterior density (not score-equivalent).

      }

    \item \emph{numerical data} (multivariate normal distribution)
      \itemize{

        \item a score-equivalent {gaussian posterior density} (\code{bge}).

      }

  }

}
\section{Whitelist and blacklist support}{

  All learning algorithms support arc whitelisting and blacklisting:

  \itemize{

    \item blacklisted arcs are never present in the graph.
    \item arcs whitelisted in one direction only (i.e.
      \eqn{A \rightarrow B}{A -> B} is whitelisted but
      \eqn{B \rightarrow A}{B -> A} is not) have the
      respective reverse arcs blacklisted, and are always
      present in the graph.
    \item arcs whitelisted in both directions (i.e. both
      \eqn{A \rightarrow B}{A -> B} and \eqn{B \rightarrow A}{B -> A}
      are whitelisted) are present in the graph,
      but their direction is set by the learning algorithm.
  }

  Any arc whitelisted and blacklisted at the same time is assumed to
  be whitelisted, and is thus removed from the blacklist.

}
\section{Error detection and correction: the strict mode}{

  Optimized implementations of the constraint-based algorithms rely
  heavily on backtracking to reduce the number of tests needed by
  the learning procedure. This approach may hide errors either
  in the Markov blanket or the neighbourhood detection phase in
  some particular cases, such as when hidden variables are present
  or there are external (logical) constraints on the interactions
  between the variables.

  On the other hand in the unoptimized implementations the Markov
  blanket and neighbour detection of each node is completely
  independent from the rest of the learning process. Thus it may
  happen that the Markov blanket or the neighbourhoods are not
  symmetric (i.e. A is in the Markov blanket of B but not vice versa),
  or that some arc directions conflict with each other.

  The \code{strict} parameter enables some measure of error
  correction, which may help to retrieve a good model
  even when the learning process would fail:

  \itemize{
    \item if \code{strict} is set to \code{TRUE}, every error
      stops the learning process and results in an error message.
    \item if \code{strict} is set to \code{FALSE}:
    \enumerate{
      \item v-structures are applied to the network structure
        in lowest-p.value order; if any arc is already oriented
        in the opposite direction, the v-structure is discarded.
      \item nodes which causes asymmetries in any Markov blanket
        are removed from that Markov blanket.
      \item nodes which causes asymmetries in any neighbourhood
        are removed from that neighbourhood.
    }
  }

}
\references{

  A. Agresti. Categorical Data Analysis. John Wiley & Sons, Inc., 2002.

  D. Margaritis. Learning Bayesian Network Model Structure from Data. PhD
    thesis, School of Computer Science, Carnegie-Mellon University, Pittsburgh,
    PA, May 2003. Available as Technical Report CMU-CS-03-153.

  J. Pearl. Probabilistic Reasoning in Intelligent Systems: Networks of Plausible
     Inference. Morgan Kaufmann Publishers Inc., 1988.

  I. Tsamardinos, C. F. Aliferis, and A. Statnikov. Algorithms for large scale Markov
     blanket discovery. In Proceedings of the Sixteenth International Florida
     Artificial Intelligence Research Society Conference, pages 376-381. AAAI Press,
     2003.

  S. Yaramakala, D. Margaritis. Speculative Markov Blanket Discovery for
     Optimal Feature Selection. In Proceedings of the Fifth IEEE International
     Conference on Data Mining, pages 809-812. IEEE Computer Society, 2005.

}
\examples{
library(bnlearn)
data(learning.test)

## Simple learning
# first try the Grow-Shrink algorithm
res = gs(learning.test)
# plot the network structure.
plot(res)
# now try the Incremental Association algorithm.
res2 = iamb(learning.test)
# plot the new network structure.
plot(res2)
# the network structures seem to be identical, don't they?
compare(res, res2)
# [1] TRUE
# how many tests each of the two algorithms used?
res$learning$ntests
# [1] 43
res2$learning$ntests
# [1] 50
# and the unoptimized implementation of these algorithms?
\dontrun{gs(learning.test, optimized = FALSE)$ntests}
# [1] 90
\dontrun{iamb(learning.test, optimized = FALSE)$ntests}
# [1] 116

## Greedy search
res = hc(learning.test)
plot(res)

## Another simple example (gaussian data)
data(gaussian.test)
# first try the Grow-Shrink algorithm
res = gs(gaussian.test)
plot(res)

## Blacklist and whitelist use
# the arc B - F should not be there?
blacklist = data.frame(from = c("B", "F"), to = c("F", "B"))
blacklist
#   from to
# 1    B  F
# 2    F  B
res3 = gs(learning.test, blacklist = blacklist)
plot(res3)
# force E - F direction (E -> F).
whitelist = data.frame(from = c("E"), to = c("F"))
whitelist
#   from to
# 1    E  F
res4 = gs(learning.test, whitelist = whitelist)
plot(res4)
# use both blacklist and whitelist.
res5 = gs(learning.test, whitelist = whitelist, blacklist = blacklist)
plot(res5)

## Debugging
# use the debugging mode to see the learning algorithms
# in action.
res = gs(learning.test, debug = TRUE)
res = hc(learning.test, debug = TRUE)
# log the learning process for future reference.
\dontrun{
sink(file = "learning-log.txt")
res = gs(learning.test, debug = TRUE)
sink()
}
# if something seems wrong, try the unoptimized version
# in strict mode (inconsistencies trigger errors):
\dontrun{
res = gs(learning.test, optimized = FALSE, strict = TRUE, debug = TRUE)
}
# or disable strict mode to try to fix errors on the fly:
\dontrun{
res = gs(learning.test, optimized = FALSE, strict = FALSE, debug = TRUE)
}

}
\keyword{ package }
