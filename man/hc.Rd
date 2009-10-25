\name{score-based algorithms}
\alias{score-based algorithms}
\alias{hc}
\title{ Score-based learning algorithms }
\description{

  Learn the structure of a Bayesian network using a
  hill-climbing (HC) greedy search.

}
\usage{
hc(x, start = NULL, whitelist = NULL, blacklist = NULL,
  score = NULL, ..., debug = FALSE, restart = 0,
  perturb = 1, max.iter = Inf, optimized = TRUE)
}
\arguments{
  \item{x}{a data frame, containing the variables in the model.}
  \item{start}{an object of class \code{bn}, the preseeded directed
      acyclic graph used to initialize the algorithm. If none is
      specified, an empty one (i.e. without any arc) is used.}
  \item{whitelist}{a data frame with two columns (optionally labeled
     "from" and "to"), containing a set of arcs to be included in the
      graph. }
  \item{blacklist}{a data frame with two columns (optionally labeled
     "from" and "to"), containing a set of arcs not to be included in
      the graph. }
  \item{score}{a character string, the label of the network score to
      be used in the algorithm. If none is specified, the default
      score is the\emph{Bayesian Information Criterion} for both
      discrete and continuous data sets. See \code{\link{bnlearn-package}}
      for details.}
  \item{\dots}{additional tuning parameters for the network score.
      See \code{\link{score}} for details.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output
       is printed; otherwise the function is completely silent.}
  \item{restart}{an integer, the number of random restarts.}
  \item{perturb}{an integer, the number of attempts to randomly
       insert/remove/reverse an arc on every random restart.}
  \item{max.iter}{an integer, the maximum number of iterations.}
  \item{optimized}{a boolean value. See \code{\link{bnlearn-package}}
       for details.}
}
\value{

  An object of class \code{bn}.
  See \code{\link{bn-class}} for details.

}
\references{

  Korb K, Nicholson AE (2003). \emph{Bayesian Artificial Intelligence}.
      Chapman & Hall/CRC.

  Margaritis D (2003). \emph{Learning Bayesian Network Model Structure
      from Data}. Ph.D. thesis, School of Computer Science, Carnegie-Mellon
      University, Pittsburgh, PA. Available as Technical Report CMU-CS-03-153.

  Daly R, Shen Q (2007). "Methods to Accelerate the Learning of Bayesian
    Network Structures". In "Proceedings of the 2007 UK Workshop on
    Computational Intelligence", Imperial College, London.

}
\author{ Marco Scutari }
\seealso{\code{\link{constraint-based algorithms}}, \code{\link{hybrid algorithms}}, \cr
  \code{\link{local discovery algorithms}}. }
\keyword{multivariate}
\keyword{models}
\keyword{graphs}
