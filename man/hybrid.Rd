\name{hybrid algorithms}
\alias{hybrid algorithms}
\alias{rsmax2}
\alias{mmhc}
\title{Hybrid structure learning algorithms}
\description{

  Learn the structure of a Bayesian network with the Max-Min Hill Climbing
  (MMHC) and the more general 2-phase Restricted Maximization (RSMAX2)
  hybrid algorithms.

}
\usage{
rsmax2(x, whitelist = NULL, blacklist = NULL, restrict = "si.hiton.pc",
  maximize = "hc", restrict.args = list(), maximize.args = list(), debug = FALSE)
mmhc(x, whitelist = NULL, blacklist = NULL, restrict.args = list(),
  maximize.args = list(), debug = FALSE)
}
\arguments{
  \item{x}{a data frame containing the variables in the model.}
  \item{whitelist}{a data frame with two columns (optionally labeled "from" and
    "to"), containing a set of arcs to be included in the graph.}
  \item{blacklist}{a data frame with two columns (optionally labeled "from" and
    "to"), containing a set of arcs not to be included in the graph.}
  \item{restrict}{a character string, the constraint-based or local search
    algorithm to be used in the \dQuote{restrict} phase. See
    \code{\link{bnlearn-package}} and the documentation of each algorithm for
    details.}
  \item{maximize}{a character string, the score-based algorithm to be used in
    the \dQuote{maximize} phase. Possible values are \code{hc} and \code{tabu}.
    See \code{\link{bnlearn-package}} for details.}
  \item{restrict.args}{a list of arguments to be passed to the algorithm
    specified by \code{restrict}, such as \code{test} or \code{alpha}.}
  \item{maximize.args}{a list of arguments to be passed to the algorithm
    specified by \code{maximize}, such as \code{restart} for hill-climbing or
    \code{tabu} for tabu search.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
}
\note{

  \code{mmhc()} is simply \code{rsmax2()} with \code{restrict} set to
  \code{mmpc} and \code{maximize} set to \code{hc}.

}
\value{

  An object of class \code{bn}. See \code{\link{bn-class}} for details.

}
\references{

  Tsamardinos I, Brown LE, Aliferis CF (2006). "The Max-Min Hill-Climbing
    Bayesian Network Structure Learning Algorithm". \emph{Machine Learning},
    \strong{65}(1):31--78.

}
\author{Marco Scutari}
\seealso{\link{local discovery algorithms},
  \link{score-based algorithms}, \link{constraint-based algorithms}.}
\keyword{structure learning}
