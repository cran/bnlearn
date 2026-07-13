\name{scm class}
\alias{scm class}
\alias{scm-class}
\alias{as.scm}
\alias{as.scm.bn}
\title{The scm class structure}
\description{

  The structure of an object of S3 class \code{scm}.

}
\usage{
as.scm(x, ...)
\method{as.scm}{bn}(x, ...)
}
\arguments{
  \item{x}{an object of class \code{bn}.}
  \item{\dots}{additional arguments, currently ignored.}
}
\details{

  An object of class \code{scm} is a list containing at least the following
  components:

  \itemize{

    \item \code{roles}: a list with three named elements:

    \itemize{

      \item \code{factual}: a vector of character strings, the labels of the
        factual nodes in the structural causal model.
      \item \code{exogenous}: a vector of character strings, the labels of the
        exogenous nodes in the structural causal model.
      \item \code{counterfactual}: a vector of character strings, the labels of
        the counterfactual nodes in the structural causal model. Unlike previous
        elements, it can be empty.

    }

    \item \code{nodes}: a list. Each element is named after a node and contains
      the following elements:

    \itemize{

      \item \code{counterfactual}: a character string.

      \itemize{
        \item For a factual node, it contains the label of the corresponding
          counterfactual node (if any).
        \item For an exogenous node, it contains the label of the counterfactual
          node linked to it (if any).
        \item For a counterfactual node, it is set to \code{character(0)}.
      }

      \item \code{exogenous}: a character string.

      \itemize{
        \item For a factual or counterfactual node, it contains the label of the
          corresponding exogenous node (if any).
        \item For an exogenous node, it is set to \code{character(0)}.
      }

      \item \code{factual}: a character string.

      \itemize{
        \item For a counterfactual node, it contains the label of the
          corresponding factual node.
        \item For an exogenous node, it contains the label of the factual node
          linked to it (if any).
        \item For a factual node, it is set to \code{character(0)}.
      }

      \item \code{parents}: the parents of the node (a vector of character
        strings). Note that exogenous nodes are not listed among the parents of
        factual and counterfactual nodes.
      \item \code{children}: the children of the node (a vector of character
        strings). Note that factual and counterfactual nodes are not listed
        among the children of exogenous nodes.

    }

    \item \code{arcs}: the arcs of the structural causal model (a two-column
      matrix, whose columns are labelled \code{from} and \code{to}), including
      the arcs between exogenous and factual/counterfactual nodes.

  }

}
\author{Marco Scutari}
\keyword{classes}
