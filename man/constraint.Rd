\name{constraint-based algorithms}
\alias{constraint-based algorithms}
\alias{gs}
\alias{iamb}
\alias{fast.iamb}
\alias{inter.iamb}
\alias{mmpc}
\alias{si.hiton.pc}
\title{Constraint-based structure learning algorithms}
\description{

  Learn the equivalence class of a directed acyclic graph (DAG) from data using
  the Grow-Shrink (GS), the Incremental Association (IAMB), the Fast Incremental
  Association (Fast-IAMB), the Interleaved Incremental Association (Inter-IAMB),
  the Max-Min Parents and Children (MMPC) or the Semi-Interleaved HITON-PC 
  constraint-based algorithms.
}
\usage{
gs(x, cluster = NULL, whitelist = NULL, blacklist = NULL, test = NULL,
  alpha = 0.05, B = NULL, debug = FALSE, optimized = TRUE, strict = FALSE,
  undirected = FALSE)
iamb(x, cluster = NULL, whitelist = NULL, blacklist = NULL, test = NULL,
  alpha = 0.05, B = NULL, debug = FALSE, optimized = TRUE, strict = FALSE,
  undirected = FALSE)
fast.iamb(x, cluster = NULL, whitelist = NULL, blacklist = NULL, test = NULL,
  alpha = 0.05, B = NULL, debug = FALSE, optimized = TRUE, strict = FALSE,
  undirected = FALSE)
inter.iamb(x, cluster = NULL, whitelist = NULL, blacklist = NULL, test = NULL,
  alpha = 0.05, B = NULL, debug = FALSE, optimized = TRUE, strict = FALSE,
  undirected = FALSE)
mmpc(x, cluster = NULL, whitelist = NULL, blacklist = NULL, test = NULL,
  alpha = 0.05, B = NULL, debug = FALSE, optimized = TRUE, strict = FALSE,
  undirected = TRUE)
si.hiton.pc(x, cluster = NULL, whitelist = NULL, blacklist = NULL, test = NULL,
  alpha = 0.05, B = NULL, debug = FALSE, optimized = TRUE, strict = FALSE,
  undirected = TRUE)
}
\arguments{
  \item{x}{a data frame containing the variables in the model.}
  \item{cluster}{an optional cluster object from package \pkg{parallel}. See
    \code{\link{parallel integration}} for details and a simple example.}
  \item{whitelist}{a data frame with two columns (optionally labeled "from"
    and "to"), containing a set of arcs to be included in the graph.}
  \item{blacklist}{a data frame with two columns (optionally labeled "from"
    and "to"), containing a set of arcs not to be included in the graph.}
  \item{test}{a character string, the label of the conditional independence
    test to be used in the algorithm. If none is specified, the default test
    statistic is the \emph{mutual information} for categorical variables, the
    Jonckheere-Terpstra test for ordered factors and the \emph{linear
    correlation} for continuous variables. See \code{\link{bnlearn-package}}
    for details.}
  \item{alpha}{a numeric value, the target nominal type I error rate.}
  \item{B}{a positive integer, the number of permutations considered for each
    permutation test. It will be ignored with a warning if the conditional
    independence test specified by the \code{test} argument is not a 
    permutation test.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
  \item{optimized}{a boolean value. See \code{\link{bnlearn-package}} for
    details.}
  \item{strict}{a boolean value. If \code{TRUE} conflicting results in the
    learning process generate an error; otherwise they result in a warning.}
  \item{undirected}{a boolean value. If \code{TRUE} no attempt will be made to
    determine the orientation of the arcs; the returned (undirected) graph
    will represent the underlying structure of the Bayesian network.}
}
\value{

  An object of class \code{bn}.
  See \code{\link{bn-class}} for details.

}
\references{

  \bold{for GS:}

  Margaritis D (2003). \emph{Learning Bayesian Network Model Structure from
      Data}. Ph.D. thesis, School of Computer Science, Carnegie-Mellon
      University, Pittsburgh, PA. Available as Technical Report CMU-CS-03-153.

  \bold{for IAMB:}

  Tsamardinos I, Aliferis CF, Statnikov A (2003). "Algorithms for Large Scale
      Markov Blanket Discovery". In "Proceedings of the Sixteenth International
      Florida Artificial Intelligence Research Society Conference", pp. 376-381.
      AAAI Press.

  \bold{for Fast-IAMB and Inter-IAMB:}

  Yaramakala S, Margaritis D (2005). "Speculative Markov Blanket Discovery for
      Optimal Feature Selection". In "ICDM '05: Proceedings of the Fifth IEEE
      International Conference on Data Mining", pp. 809-812. IEEE Computer
      Society.

  \bold{for MMPC:}

  Tsamardinos I, Aliferis CF, Statnikov A (2003). "Time and Sample
      Efficient Discovery of Markov Blankets and Direct Causal Relations".
      In "KDD '03: Proceedings of the Ninth ACM SIGKDD International
      Conference on Knowledge Discovery and Data Mining", pp. 673-678. ACM.

  Tsamardinos I, Brown LE, Aliferis CF (2006). "The Max-Min Hill-Climbing
      Bayesian Network Structure Learning Algorithm". \emph{Machine Learning},
      \strong{65}(1), 31-78.

  \bold{for the Semi-Interleaved HITON-PC:}

  Aliferis FC, Statnikov A, Tsamardinos I, Subramani M, Koutsoukos XD (2010).
      "Local Causal and Markov Blanket Induction for Causal Discovery and
      Feature Selection for Classification Part I: Algorithms and Empirical
      Evaluation". \emph{Journal of Machine Learning Research}, \strong{11},
      171-234.

}
\author{ Marco Scutari }
\seealso{\link{local discovery algorithms},
  \link{score-based algorithms}, \link{hybrid algorithms}.}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}
