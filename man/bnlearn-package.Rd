\name{bnlearn-package}
\alias{bnlearn-package}
\alias{bnlearn}
\docType{package}
\title{

  Bayesian network structure learning.

}
\description{

  Bayesian network learning via constraint-based, score-based and hybrid
  algorithms.

}
\details{

  \tabular{ll}{

    Package: \tab bnlearn\cr
    Type: \tab Package\cr
    Version: \tab 1.8\cr
    Date: \tab 2010-01-05\cr
    License: \tab GPLv2 or later

  }

  This package implements some algorithms for learning
  the structure of Bayesian networks.

  \emph{Constraint-based algorithms}, also known as \emph{conditional
  independence learners}, are all optimized derivatives of the
  \emph{Inductive Causation} algorithm (Verma and Pearl, 1991).
  These algorithms use conditional independence tests to detect the Markov
  blankets of the variables, which in turn are used to compute the
  structure of the Bayesian network.

  \emph{Score-based learning algorithms} are general purpose heuristic
  optimization algorithms which rank network structures with respect to
  a goodness-of-fit score.

  \emph{Hybrid algorithms} combine aspects of both constraint-based and
   score-based algorithms, as they use conditional independence tests
   (usually to reduce the search space) and network scores (to find the
   optimal network in the reduced space) at the same time.

}
\author{

  Marco Scutari\cr
  Department of Statistical Sciences\cr
  University of Padova

  Maintainer: Marco Scutari \email{marco.scutari@gmail.com}

}
\section{Available constraint-based learning algorithms}{

  \itemize{

    \item \emph{Grow-Shrink} (\code{\link{gs}}): based on the \emph{Grow-Shrink
      Markov Blanket}, the first (and simplest) Markov blanket detection algorithm
      (Margaritis, 2003) used in a structure learning algorithm.
    \item \emph{Incremental Association} (\code{\link{iamb}}): based on the
      Markov blanket detection algorithm of the same name (Tsamardinos et
      al., 2003), which is based on a two-phase selection scheme (a forward
      selection followed by an attempt to remove false positives).
    \item \emph{Fast Incremental Association} (\code{\link{fast.iamb}}): a
      variant of IAMB which uses speculative stepwise forward selection to
      reduce the number of conditional independence tests (Yaramakala and
      Margaritis, 2005).
    \item \emph{Interleaved Incremental Association} (\code{\link{inter.iamb}}):
      another variant of IAMB which uses forward stepwise selection (Tsamardinos
      et al., 2003) to avoid false positives in the Markov blanket detection
      phase.

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

  The computational complexity of these algorithms is polynomial in the 
  number of tests, usually \eqn{O(N^2)}{O(N^2)} (\eqn{O(N^4)}{O(N^4)}
  in the worst case scenario), where \eqn{N}{N} is the number of variables. 
  Execution time scales linearly with the size of the data set.

}
\section{Available score-based learning algorithms}{

  \itemize{

    \item \emph{Hill-Climbing} (\code{\link{hc}}): a \emph{hill climbing}
      greedy search on the space of the directed graphs. The optimized
      implementation uses score caching, score decomposability and score
      equivalence to reduce the number of duplicated tests.
    \item \emph{Tabu Search} (\code{\link{tabu}}): a modified hill
      climbing able to escape local optima by selecting a network that
      minimally decreases the score function.

  }

  Random restart with a configurable number of perturbing operations 
  is implemented for both algorithms.

}
\section{Available hybrid learning algorithms}{

  \itemize{

    \item \emph{Max-Min Hill-Climbing} (\code{\link{mmhc}}): a hybrid
      algorithm which combines the Max-Min Parents and Children algorithm
      (to restrict the search space) and the Hill-Climbing algorithm
      (to find the optimal network structure in the restricted space).
    \item \emph{Restricted Maximization} (\code{\link{rsmax2}}): a more
      general implementation of the Max-Min Hill-Climbing, which can use
      any combination of constraint-based and score-based algorithms.

  }

}
\section{Other (constraint-based) local discovery algorithms}{

  These algorithms learn the structure of the undirected graph underlying
  the Bayesian network, which is known as the \emph{skeleton} of the
  network or the \emph{correlation graph}. Therefore all the arcs are
  undirected, and no attempt is made to detect their orientation. They
  are often used in hybrid learning algorithms.

  \itemize{

    \item \emph{Max-Min Parents and Children} (\code{\link{mmpc}}): a
      forward selection technique for neighbourhood detection based on
      the maximization of the minimum association measure observed with
      any subset of the nodes selected in the previous iterations
      (Tsamardinos, Brown and Aliferis, 2006).

  }

  All these algorithms have three implementations (unoptimized,
  optimized and cluster-aware) like other constraint-based algorithms.

}
\section{Available (conditional) independence tests}{

  The conditional independence tests used in \emph{constraint-based}
  algorithms in practice are statistical tests on the data set. Available
  tests (and the respective labels) are:

  \itemize{

    \item \emph{discrete case} (multinomial distribution)
      \itemize{

        \item \emph{mutual information}: an information-theoretic
          distance measure. It's proportional to the log-likelihood ratio
          (they differ by a \eqn{2n}{2n} factor) and is related to the
          deviance of the tested models. Both the asymptotic \eqn{\chi^2}{chi-square}
          test (\code{mi}) and the Monte Carlo permutation test (\code{mc-mi})
          are implemented.
        \item \emph{shrinkage estimator} for the \emph{mutual information}
          (\code{mi-sh}): an improved asymptotic \eqn{\chi^2}{chi-square} test
          based on the James-Stein estimator for the mutual information.
        \item \emph{Pearson's \eqn{X^2}{X^2}}: the classical Pearson's
          \eqn{X^2}{X^2} test for contingency tables. Both the asymptotic
          \eqn{\chi^2}{chi-square} test (\code{x2}) and the Monte Carlo
          permutation test (\code{mc-x2}) are implemented.
        \item \emph{Akaike Information Criterion} (\code{aict}): an experimental
          AIC-based independence test, computed comparing the mutual information
          and the expected information gain.

      }
    \item \emph{continuous case} (multivariate normal distribution)
      \itemize{

        \item \emph{linear correlation}: linear correlation. Both the exact
          Student's t test (\code{cor}) and the Monte Carlo permutation test
          (\code{mc-cor}) are implemented.
        \item \emph{Fisher's Z}: a transformation of the linear correlation
          with asymptotic normal distribution. Used by commercial software
          (such as TETRAD II) for the PC algorithm (an R implementation is
          present in the \code{pcalg} package on CRAN). Both the asymptotic
          normal (\code{zf}) and the Monte Carlo permutation test (\code{mc-zf})
          are implemented.
        \item \emph{mutual information}: an information-theoretic distance
          measure. Again it's proportional to the log-likelihood ratio (they
          differ by a \eqn{2n}{2n} factor). Both the asymptotic \eqn{\chi^2}{chi-square}
          test (\code{mi-g}) and the Monte Carlo permutation test (\code{mc-mi-g})
          are implemented.

      }

  }

}
\section{Available network scores}{

  Available scores (and the respective labels) are:

  \itemize{

    \item \emph{discrete case} (multinomial distribution)
      \itemize{

        \item the multinomial \emph{log-likelihood} (\code{loglik}) score,
          which is equivalent to the \emph{entropy measure} used in Weka.
        \item the \emph{Akaike Information Criterion} score (\code{aic}).
        \item the \emph{Bayesian Information Criterion} score (\code{bic}),
          which is equivalent to the \emph{Minimum Description Length} (MDL)
          and is also known as \emph{Schwarz Information Criterion}.
        \item the logarithm of the \emph{Bayesian Dirichlet equivalent}
          score (\code{bde}), a score equivalent Dirichlet posterior density.
        \item the logarithm of the \emph{K2} score (\code{k2}), a Dirichlet
          posterior density (not score equivalent).

      }

    \item \emph{continuous case} (multivariate normal distribution)
      \itemize{

        \item the multivariate Gaussian \emph{log-likelihood} (\code{loglik-g})
          score.
        \item the corresponding \emph{Akaike Information Criterion} score (\code{aic-g}).
        \item the corresponding \emph{Bayesian Information Criterion} score (\code{bic-g}).
        \item a score equivalent \emph{Gaussian posterior density} (\code{bge}).
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

  Optimized implementations of constraint-based algorithms rely
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
  even when the learning process would otherwise fail:

  \itemize{
    \item if \code{strict} is set to \code{TRUE}, every error
      stops the learning process and results in an error message.
    \item if \code{strict} is set to \code{FALSE}:
    \enumerate{
      \item v-structures are applied to the network structure
        in lowest-p.value order; if any arc is already oriented
        in the opposite direction, the v-structure is discarded.
      \item nodes which cause asymmetries in any Markov blanket
        are removed from that Markov blanket; they are treated
        as false positives.
      \item nodes which cause asymmetries in any neighbourhood
        are removed from that neighbourhood; again they are treated
        as false positives (see Tsamardinos, Brown and Aliferis, 2006).
    }
  }

}
\references{

  (a BibTeX file with all the references cited throughout this manual
  is present in the \file{bibtex} directory of this package)

  Agresti A (2002). \emph{Categorical Data Analysis}. Wiley Series in
      Probability and Statistics. Wiley-Interscience, 2nd edition.

  Korb K, Nicholson AE (2003). \emph{Bayesian Artificial Intelligence}.
      Chapman & Hall/CRC.

  Margaritis D (2003). \emph{Learning Bayesian Network Model Structure
      from Data}. Ph.D. thesis, School of Computer Science, Carnegie-Mellon
      University, Pittsburgh, PA. Available as Technical Report CMU-CS-03-153.

  Pearl J (1988). \emph{Probabilistic reasoning in intelligent systems:
      networks of plausible inference}. Morgan Kaufmann.

  Tsamardinos I, Aliferis CF, Statnikov A (2003). "Algorithms for Large
      Scale Markov Blanket Discovery". In "Proceedings of the Sixteenth
      International Florida Artificial Intelligence Research Society
      Conference", pp. 376-381. AAAI Press.

  Tsamardinos I, Brown LE, Aliferis CF (2006). "The Max-Min Hill-Climbing
      Bayesian Network Structure Learning Algorithm". Machine Learning,
      65(1), 31-78.

  Yaramakala S, Margaritis D (2005). "Speculative Markov Blanket Discovery
      for Optimal Feature Selection". In "ICDM '05: Proceedings of the Fifth
      IEEE International Conference on Data Mining", pp. 809-812. IEEE
      Computer Society.

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
\dontrun{gs(learning.test, optimized = FALSE)$learning$ntests}
# [1] 93
\dontrun{iamb(learning.test, optimized = FALSE)$learning$ntests}
# [1] 116

## Greedy search
res = hc(learning.test)
plot(res)

## Another simple example (Gaussian data)
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
# or disable strict mode to let the algorithm fix errors on the fly:
\dontrun{
res = gs(learning.test, optimized = FALSE, strict = FALSE, debug = TRUE)
}

}
\keyword{ package }
