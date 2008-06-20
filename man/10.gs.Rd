\name{gs}
\alias{gs}
\title{ Grow-Shrink (GS) learning algorithm }
\description{

  Estimate the equivalence class of a directed acyclic graph (DAG)
  from data using the Grow-Shrink (GS) constraint-based algorithm.

}
\usage{
  gs(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
    test = NULL, alpha = 0.05, debug = FALSE, optimized = TRUE,
    strict = TRUE, direction = FALSE)
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
      independence test to be used in the algorithm. Possible
      values are \code{mi} (\emph{mutual information}), \code{mh}
      (\emph{Cochran-Mantel-Haenszel}), \code{fmi} (\emph{fast
      mutual information}), \code{aict} (\emph{AIC}), \code{cor}
      (\emph{linear correlation}), \code{zf} (\emph{Fisher's Z}).
      If none is specified, the default one is \emph{mutual
      information} for discrete data sets and \emph{linear correlation}
      for continuous ones. See \code{\link{bnlearn-package}} for details.}
  \item{alpha}{a numerical value, the target nominal type I error rate.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output
       is printed; otherwise the function is completely silent.}
  \item{optimized}{a boolean value. See \code{\link{bnlearn-package}}
       for details.}
  \item{strict}{a boolean value. If \code{TRUE} conflicting results in
      the learning process generate an error; otherwise they result
      in a warning.}
  \item{direction}{a boolean value. If \code{TRUE} each possible direction
      of each undirected arc is tested, and the one with the lowest
      p-value is accepted as the true direction for that arc.}
}
\value{

  An object of class \code{bn}.
  See \code{\link{bn-class}} for details.

}
\references{

  D. Margaritis. Learning Bayesian Network Model Structure from Data. PhD
    thesis, School of Computer Science, Carnegie-Mellon University, Pittsburgh,
    PA, May 2003. Available as Technical Report CMU-CS-03-153.

}
\author{ Marco Scutari }
\seealso{\code{\link{iamb}}, \code{\link{fast.iamb}}, \code{\link{inter.iamb}},
  \code{\link{hc}}. }
\keyword{multivariate}
\keyword{models}
\keyword{graphs}
