\name{bnlearn-package}
\alias{bnlearn-package}
\alias{bnlearn}
\docType{package}
\title{Bayesian network structure learning, parameter learning and inference}
\description{

  Bayesian network structure learning (via constraint-based, score-based and
  hybrid algorithms), parameter learning (via ML and Bayesian estimators) and
  inference (via approximate inference algorithms).

}
\details{

  \pkg{bnlearn} implements key algorithms covering all stages of Bayesian
  network modelling: data preprocessing, structure learning combining data and
  expert/prior knowledge, parameter learning, and inference (including causal
  inference via do-calculus). \pkg{bnlearn} aims to be a one-stop shop for
  Bayesian networks in R, providing the tools needed for learning and working
  with discrete Bayesian networks, Gaussian Bayesian networks and conditional
  linear Gaussian Bayesian networks on real-world data. Incomplete data with
  missing values are also supported. Furthermore the modular nature of
  \pkg{bnlearn} makes it easy to use it for simulation studies.

  Implemented structure learning algorithms include:
  \itemize{

    \item \emph{Constraint-based algorithms}, which use conditional independence
      tests to learn conditional independence constraints from data. The
      constraints in turn are used to learn the structure of the Bayesian
      network under the assumption that conditional independence implies
      graphical separation (so, two variables that are independent cannot be
      connected by an arc).


    \item \emph{Score-based algorithms}, which are general-purpose optimization
      algorithms that rank network structures with respect to a goodness-of-fit
      score.

    \item \emph{Hybrid algorithms} combine aspects of both constraint-based and
	  score-based algorithms, as they use conditional independence tests
	  (usually to reduce the search space) and network scores (to find the
	  optimal network in the reduced space) at the same time.

  }

  For more details about structure learning algorithms see
  \link{structure learning}; available conditional independence tests are
  described in \link{independence tests} and available network scores are
  described in \link{network scores}. Specialized algorithms to learn the
  structure of Bayesian network classifiers are described in
  \link{network classifiers}. All algorithms support the use of whitelists and
  blacklists to include and exclude arcs from the networks (see
  \link{whitelists and blacklists}); and many have parallel implementation
  built on the \pkg{parallel} package. Bayesian network scores support the use
  of graphical priors.

  Parameter learning approaches include both frequentist and Bayesian
  estimators. Inference is implemented using approximate algorithms via particle
  filters approaches such as likelihood weighting, and covers conditional
  probability queries, prediction and imputation.

  Additional facilities include support for bootstrap and cross-validation;
  advanced plotting capabilities implemented on top of \pkg{Rgraphviz} and
  \pkg{lattice}; model averaging; random graphs and random samples generation;
  import/export functions to integrate \pkg{bnlearn} with software such as
  Hugin and GeNIe; an associated Bayesian network repository of golden-standard
  networks at \url{https://www.bnlearn.com/bnrepository/}.

  Use \code{citation("bnlearn")} to find out how to cite \pkg{bnlearn} in
  publications and other materials; and visit \url{https://www.bnlearn.com/} for
  more examples and code from publications using \pkg{bnlearn}.

}
\author{

  Marco Scutari\cr
  Istituto Dalle Molle di Studi sull'Intelligenza Artificiale (IDSIA)\cr

  Maintainer: Marco Scutari \email{marco.scutari@gmail.com}

}
\references{

  \bold{reference books:}

  Koller D, Friedman N (2009). \emph{Probabilistic Graphical Models: Principles
    and Techniques}. MIT Press.

  Korb K, Nicholson AE (2010). \emph{Bayesian Artificial Intelligence}.
    Chapman & Hall/CRC, 2nd edition.

  Pearl J (1988). \emph{Probabilistic Reasoning in Intelligent Systems:
    Networks of Plausible Inference}. Morgan Kaufmann.

  \bold{from the author:}

  Nagarajan R, Scutari M, Lebre S (2013). "Bayesian Networks in R with
    Applications in Systems Biology". Springer.

  Scutari M (2010). "Learning Bayesian Networks with the bnlearn R Package".
    \emph{Journal of Statistical Software}, \strong{35}(3):1--22.

  Scutari M (20107). "Bayesian Network Constraint-Based Structure Learning
    Algorithms: Parallel and Optimized Implementations in the bnlearn R
    Package". \emph{Journal of Statistical Software}, \strong{77}(2):1--20.

}
\examples{
## the workflow of Bayesian network modelling in bnlearn:
# choose the data set to work on...
data(learning.test)
# ... choose an algorithm and learn the structure of the network from the data...
net = hc(learning.test)
# ... plot it...
\dontrun{graphviz.plot(net)}
# ... learn the parameters of the network...
bn = bn.fit(net, learning.test)
# ... explore the network with a classic barchart...
\dontrun{graphviz.chart(bn)}
# ... and perform inference to answer any question that interests you!
cpquery(bn, event = (A == "a"), evidence = (C == "a"))
}
\keyword{package}
