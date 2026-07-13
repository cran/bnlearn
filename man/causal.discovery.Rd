\name{causal discovery algorithms}
\alias{causal discovery algorithms}
\alias{direct.lingam}
\title{Causal discovery algorithms}
\description{

  Learn the structure of a (causal) Bayesian network with DirectLiNGAM.

}
\usage{
direct.lingam(x, cluster, whitelist, blacklist, mi, maximize = "alasso",
  maximize.args = list(), debug = FALSE)
}
\arguments{
  \item{x}{a data frame containing the variables in the model.}
  \item{cluster}{an optional cluster object from package \pkg{parallel}.}
  \item{whitelist}{a data frame with two columns (optionally labelled "from" and
    "to"), containing a set of arcs to be included in the graph.}
  \item{blacklist}{a data frame with two columns (optionally labelled "from" and
    "to"), containing a set of arcs not to be included in the graph.}
  \item{mi}{a character string, the label of the mutual information measure used
    to identify the topological ordering of the nodes. It can take values
    \code{"pwling"} (the default) or \code{"gkernel"}.}
  \item{maximize}{a character value, the label of the method used to identify
    the parents of each node given the causal ordering. The only possible
    value is "alasso" (adaptive LASSO).}
  \item{maximize.args}{a list of arguments to be passed to the method specified
    by \code{maximize}. See below for details.}
  \item{debug}{a boolean value. If \code{TRUE}, a lot of debugging output is
    printed. Otherwise, the function is completely silent.}
}
\details{

  The \code{"alasso"} method has the following optional arguments:
  \itemize{
    \item \code{gamma}: the coefficient of the weights used in the adaptive
      LASSO. The default value is \code{1}.
    \item \code{lambda.min.ratio} and \code{pmax}: identical to the arguments
      of the same names in \code{glmnet()}, with the same default values and
      interpretation.
    \item \code{k}: the coefficient of the penalised likelihood score used to
      choose the optimal shrinkage in both ridge and LASSO models in the
      adaptive LASSO. It defaults to the BIC penalty.
  }

}
\note{

  See \code{\link{structure learning}} for a complete list of structure learning
  algorithms with the respective references.

}
\value{

  An object of class \code{bn}. See \code{\link{bn-class}} for details.

}
\author{Marco Scutari}
\seealso{\link{local discovery algorithms},
  \link{score-based algorithms}, \link{constraint-based algorithms},
  \link{hybrid algorithms}.}
\keyword{structure learning}
\keyword{causality}
