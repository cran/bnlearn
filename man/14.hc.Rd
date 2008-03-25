\name{hc}
\alias{hc}
\title{ Hill-Climbing (HC) learning algorithm }
\description{

  Learn the structure of a Bayesian network using a
  hill-climbing (HC) greedy search.

}
\usage{
hc(x, start = NULL, whitelist = NULL, blacklist = NULL,
    score = NULL, ..., debug = FALSE, optimized = TRUE)
}
\arguments{
  \item{x}{a data frame, containing the variables in the model.}
  \item{start}{an object of class \code{bn}, the preseeded network
      used to initialize the algorithm. If none is specified, an
      empty one (i.e. without any arc) is used.}
  \item{whitelist}{a data frame with two columns (optionally labeled
     "from" and "to"), containing a set of arcs to be included in the
      graph. }
  \item{blacklist}{a data frame with two columns (optionally labeled
     "from" and "to"), containing a set of arcs not to be included in
      the graph. }
  \item{score}{a character string, the label of the network score to
      be used in the algorithm. Possible values are \code{lik}
      (\emph{multinomial likelihood}), \code{loglik} (\emph{multinomial
      loglikelihood}), \code{aic} (\emph{Akaike Information Criterion}),
      \code{bic} (\emph{Bayesian Information Criterion}), \code{k2}
      (\emph{K2 score}), \code{bde} or \code{dir} (\emph{Bayesian
      dirichlet posterior density}), \code{bge} (\emph{Bayesian
      gaussian posterior density}). If none is specified, the default
      score is the \emph{Akaike Information Criterion} for discrete
      data sets and the \emph{Bayesian gaussian posterior density}
      for discrete ones. See \code{\link{bnlearn-package}}
      for details.}
  \item{...}{additional tuning parameters for the network score.
      See \code{\link{score}} for details.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output
       is printed; otherwise the function is completely silent.}
  \item{optimized}{a boolean value. See \code{\link{bnlearn-package}}
       for details.}
}
\value{

  An object of class \code{bn}.
  See \code{\link{bn-class}} for details.

}
\references{

  K. Korb and A. Nicholson. Bayesian artificial intelligence. Chapman and
   Hall, 2004.

  D. Margaritis. Learning Bayesian Network Model Structure from Data. PhD
    thesis, School of Computer Science, Carnegie-Mellon University, Pittsburgh,
    PA, May 2003. Available as Technical Report CMU-CS-03-153.

  R. Daly and Q. Shen. Methods to accelerate the learning of bayesian network
    structures. Proceedings of the 2007 UK Workshop on Computational
    Intelligence.

}
\author{ Marco Scutari }
\seealso{\code{\link{gs}}, \code{\link{iamb}}, \code{\link{fast.iamb}},
  \code{\link{inter.iamb}}. }
\keyword{multivariate}
\keyword{models}
\keyword{graphs}
