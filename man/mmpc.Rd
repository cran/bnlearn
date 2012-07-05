\name{local discovery algorithms}
\alias{local discovery algorithms}
\alias{mmpc}
\alias{si.hiton.pc}
\alias{aracne}
\alias{chow.liu}
\title{ Local discovery structure learning algorithms }
\description{

  Learn the skeleton of a directed acyclic graph (DAG) from data using
  the Max-Min Parents and Children (MMPC) and the Semi-Interleaved
  HITON-PC constraint-based algorithms. ARACNE and Chow-Liu learn an 
  approximation of that structure using pairwise mutual information
  coefficients.

}
\usage{
mmpc(x, cluster = NULL, whitelist = NULL, blacklist = NULL, test = NULL,
  alpha = 0.05, B = NULL, debug = FALSE, optimized = TRUE, strict = FALSE)
si.hiton.pc(x, cluster = NULL, whitelist = NULL, blacklist = NULL, test = NULL,
  alpha = 0.05, B = NULL, debug = FALSE, optimized = TRUE, strict = FALSE)

aracne(x, whitelist = NULL, blacklist = NULL, mi = NULL, debug = FALSE)
chow.liu(x, whitelist = NULL, blacklist = NULL, mi = NULL, debug = FALSE)
}
\arguments{
  \item{x}{a data frame containing the variables in the model.}
  \item{cluster}{an optional cluster object from package \pkg{snow}.
      See \code{\link{snow integration}} for details and a simple
      example.}
  \item{whitelist}{a data frame with two columns (optionally labeled
     "from" and "to"), containing a set of arcs to be included in the
      graph. }
  \item{blacklist}{a data frame with two columns (optionally labeled
     "from" and "to"), containing a set of arcs not to be included in
      the graph. }
  \item{mi}{a character string, the estimator used for the pairwise
      (i.e. unconditional) mutual information coefficients in the ARACNE
      and Chow-Liu algorithms. Possible values are \code{mi} (discrete
      mutual information) and \code{mi-g} (Gaussian mutual information).}
  \item{test}{a character string, the label of the conditional
      independence test to be used in the algorithm. If none is
      specified, the default test statistic is the \emph{mutual information}
      for discrete data sets and the \emph{linear correlation} for
      continuous ones. See \code{\link{bnlearn-package}} for details.}
  \item{alpha}{a numeric value, the target nominal type I error rate.}
  \item{B}{a positive integer, the number of permutations considered
      for each permutation test. It will be ignored with a warning if
      the conditional independence test specified by the \code{test}
      argument is not a permutation test.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output
      is printed; otherwise the function is completely silent.}
  \item{optimized}{a boolean value. See \code{\link{bnlearn-package}}
      for details.}
  \item{strict}{a boolean value. If \code{TRUE} conflicting results in
      the learning process generate an error; otherwise they result in
      a warning.}
}
\value{

  An object of class \code{bn}.
  See \code{\link{bn-class}} for details.

}
\references{

  Tsamardinos I, Aliferis CF, Statnikov A (2003). "Time and Sample
      Efficient Discovery of Markov Blankets and Direct Causal Relations".
      In "KDD '03: Proceedings of the Ninth ACM SIGKDD International
      Conference on Knowledge Discovery and Data Mining", pp. 673-678. ACM.

  Tsamardinos I, Brown LE, Aliferis CF (2006). "The Max-Min Hill-Climbing
      Bayesian Network Structure Learning Algorithm". \emph{Machine Learning},
      \strong{65}(1), 31-78.

  Aliferis FC, Statnikov A, Tsamardinos I, Subramani M, Koutsoukos XD (2010).
      "Local Causal and Markov Blanket Induction for Causal Discovery and
      Feature Selection for Classification Part I: Algorithms and Empirical
      Evaluation". \emph{Journal of Machine Learning Research}, \strong{11},
      171-234.

  Margolin AA, Nemenman I, Basso K, Wiggins C, Stolovitzky G, Dalla Favera R,
      Califano A (2006). "ARACNE: An Algorithm for the Reconstruction of
      Gene Regulatory Networks in a Mammalian Cellular Context". \emph{BMC
      Bioinformatics}, \strong{7}(Suppl 1):S7.

}
\author{ Marco Scutari }
\seealso{\link{constraint-based algorithms}, \link{score-based algorithms},
  \link{hybrid algorithms}. }
\keyword{multivariate}
\keyword{models}
\keyword{graphs}
