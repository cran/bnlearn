\name{mmpc}
\alias{mmpc}
\title{ Max-Min Parents and Children (MMPC) learning algorithm }
\description{

  Estimate the underlying structure of a directed acyclic graph (DAG)
  from data using the Max-Min Parents and Children (MMPC)
  constraint-based algorithm.

}
\usage{
  mmpc(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
    test = NULL, alpha = 0.05, debug = FALSE, optimized = TRUE,
    strict = FALSE)
}
\arguments{
  \item{x}{a data frame, containing the variables in the model.}
  \item{cluster}{an optional cluster object from package \pkg{snow}.
      See \code{\link{snow integration}} for details and a simple
      example.}
  \item{whitelist}{a data frame with two columns (optionally labeled
     "from" and "to"), containing a set of arcs to be included in the
      graph. }
  \item{blacklist}{a data frame with two columns (optionally labeled
     "from" and "to"), containing a set of arcs not to be included in
      the graph. }
  \item{test}{a character string, the label of the conditional
      independence test to be used in the algorithm. If none is
      specified, the default test statistic is the \emph{mutual information}
      for discrete data sets and the \emph{linear correlation} for 
      continuous ones. See \code{\link{bnlearn-package}} for details.}
  \item{alpha}{a numeric value, the target nominal type I error rate.}
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

  I. Tsamardinos, C. F. Aliferis, A. Statnikov. Time and sample
     efficient discovery of Markov blankets and direct causal
     relations. Proceedings of the Ninth International Conference
     on Knowledge Discovery and Data Mining KDD, pages 673-8, 2003.

  I. Tsamardinos, L. E. Brown, C. Aliferis. The max-min hill-climbing
     Bayesian network learning algorithm. Machine Learning, 65(1),
     pages 31-78. Kluwer Academic Publishers, 2006.

}
\author{ Marco Scutari }
\seealso{\code{\link{gs}}, \code{\link{fast.iamb}}, \code{\link{iamb}},
  \code{\link{inter.iamb}}, \code{\link{hc}}. }
\keyword{multivariate}
\keyword{models}
\keyword{graphs}
