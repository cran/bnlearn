\name{causal inference}
\alias{mutilated}
\alias{twin}
\alias{intervention}
\alias{counterfactual}
\title{Perform causal inference}
\description{

  Perform interventions and counterfactuals, using Bayesian networks as causal
  models.

}
\usage{
# network transformations used in causal inference.
mutilated(x, evidence)
twin(x)

# causal inference types.
intervention(x, evidence)
counterfactual(x, evidence, merging = TRUE)
}
\arguments{
  \item{x}{an object of class \code{bn}, \code{bn.fit} or \code{scm}.}
  \item{evidence}{a list, in the same format as the evidence for likelihood
    weighting (see \code{\link{cpquery}}).}
  \item{merging}{a boolean value. If \code{TRUE}, the twin network is reduced
    by node merging after introducing the counterfactual.}
}
\details{

  \code{mutilated()} constructs the mutilated network arising from an ideal
  intervention setting the nodes involved to the values specified by
  \code{evidence}.

  \code{twin()} constructs the twin network required to compute counterfactuals,
  separating noise (exogenous) from functional (factual) nodes and adding
  counterfactual nodes.

  \code{intervention()} is an alias of \code{mutilated()}.

  \code{counterfactual()} constructs the twin network and introduces a
  counterfactual (that is, an intervention on the counterfactual nodes) before
  applying node merging (optional). The intervention in \code{evidence} is
  specified using the labels of the nodes in \code{x}, and is then applied to
  the corresponding counterfactual nodes.

  Both \code{twin()} and \code{counterfactual()} label nodes as follows:
  \itemize{

    \item A counterfactual node has the label of the corresponding factual node
      with a period appended (as in \code{paste0(label, ".")}).
    \item A noise node has the label of the corresponding factual node with a
      \code{"u"} prepended (as in \code{paste0("u", label)}).

  }

  Evidence in \code{intervention()} can reference factual and counterfactual
  nodes. Evidence in \code{counterfactual()} can only reference counterfactual
  nodes; factual nodes will be replaced by the corresponding counterfactual
  nodes with a warning. Referencing exogenous nodes in either function will
  produce an error.

}
\value{

  A \code{bn}, \code{bn.fit} or \code{scm} object, depending on the class of
  \code{x}.

  Objects returned by \code{twin()} will additionally have class \code{bn.twin}
  (if \code{x} is a \code{bn} or \code{bn.fit} object) or \code{scm.twin} (if
  \code{x} is an \code{scm} object).

  Objects returned by \code{counterfactual()} will additionally have classes
  \code{bn.twin} and \code{bn.ctf} (if \code{x} is a \code{bn} or \code{bn.fit}
  object) or \code{scm.twin} (if \code{x} is an \code{scm} object).

}
\references{

  Balke A, Pearl J (2022). "Probabilistic Evaluation of Counterfactual Queries."
    \emph{Probabilistic and Causal Inference: The Works of Judea Pearl},
    237--254. The original is in the Proceedings of AAAI, 1994.

  Graham D, Lee CM, Perov Y (2019). "Copy, Paste, Infer: A Robust Analysis of
    Twin Networks for Counterfactual Inference." \emph{NeurIPS19 CausalML
    workshop}.

  Koller D, Friedman N (2009). \emph{Probabilistic Graphical Models: Principles
    and Techniques}. MIT Press.

  Shpitser I, PearlJ (2007). "What counterfactuals can be tested."
    \emph{Proceedings of the 23rd Conference on Uncertainty in Artificial
    Intelligence (UAI)}, 352--359.

}
\examples{
# example from Graham, Lee and Perov.
dag = model2network("[C][A|C][B|C][D|A:B]")
intervention(dag, evidence = list(A = 0))
twin(dag)
counterfactual(dag, evidence = list(A. = 0), merging = FALSE)
counterfactual(dag, evidence = list(A. = 0), merging = TRUE)
}
\author{Marco Scutari}
\keyword{inference}
\keyword{causality}
