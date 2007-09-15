\name{fast.iamb}
\alias{fast.iamb}
\title{ Fast Incremental Association (Fast-IAMB) learning algorithm }
\description{

  Estimate the equivalence class of a directed acyclic graph (DAG)
  from data using the Fast Incremental Association (FAST-IAMB)
  Constraint-based algorithm.

}
\usage{

  fast.iamb(x, cluster = NULL, whitelist = NULL, blacklist = NULL, 
    test = "mi", alpha = 0.05, debug = FALSE, optimized = TRUE, 
    strict = TRUE, direction = FALSE)

}
\arguments{
  \item{x}{a data frame, containing the variables in the model.}
  \item{cluster}{See \code{\link{bnlearn-package}} for details.}
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
      mutual information}), \code{cor} (\emph{linear correlation}),
      \code{zf} (\emph{Fisher's Z}). See
      \code{\link{bnlearn-package}} for details.}
  \item{alpha}{a numerical value, the target nominal type I error rate.}
  \item{debug}{a boolean value. If TRUE a lot of debugging output is
       printed; otherwise the function is completely silent.}
  \item{optimized}{a boolean value. See \code{\link{bnlearn-package}} 
       for details.}
  \item{strict}{a boolean value. If TRUE conflicting results in the
      structure learning generate an error; otherwise they result in
      a warning.}
  \item{direction}{a boolean value. If TRUE each possible direction 
      of each undirected arc is tested, and the one with the lowest
      p-value is accepted as the true direction for that arc.}
}
\value{

  An object of class \code{bn}.
  See \code{\link{bnlearn-package}} for details.

}
\references{ 

  S. Yaramakala, D. Margaritis. Speculative Markov Blanket Discovery for
     Optimal Feature Selection. In Proceedings of the Fifth IEEE International 
     Conference on Data Mining, pages 809-812. IEEE Computer Society, 2005.

}
\author{ Marco Scutari }
\seealso{\code{\link{gs}}, \code{\link{fast.iamb}}, \code{\link{iamb}}. }
\keyword{multivariate}
\keyword{models}
\keyword{graphs}
