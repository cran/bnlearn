\name{network-classifiers}
\alias{network-classifiers}
\alias{network classifiers}
\title{Bayesian network Classifiers}
\description{

  Structure learning algorithms for Bayesian network classifiers.

}
\details{

  The algorithms are aimed at classification, and favour predictive power over
  the ability to recover the correct network structure. The implementation in
  \pkg{bnlearn} assumes that all variables, including the classifiers, are
  discrete.

  \itemize{

    \item \emph{Naive Bayes} (\code{\link{naive.bayes}}): a very simple
      algorithm that assumes all classifiers are independent and uses the target
      variable's posterior probability for classification.
    \item \emph{Tree-Augmented Naive Bayes} (\code{\link{tree.bayes}}):
      improves over Naive Bayes by using the Chow-Liu algorithm to approximate
      the dependence structure of the classifiers.

     Friedman N, Geiger D, Goldszmit M (1997). "Bayesian Network Classifiers."
       \emph{Machine Learning}, \strong{29}:131--163.

  }

}
\keyword{package}
