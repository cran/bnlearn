\name{bn.fit class}
\alias{bn.fit class}
\alias{bn.fit-class}
\alias{bn.fit.dnode}
\alias{bn.fit.gnode}
\title{

  The bn.fit class structure

}
\description{

  The structure of an object of the \code{bn.fit} S3 class.

}
\details{

  An object of class \code{bn.fit} is a list whose elements
  correspond to the nodes of the Bayesian network. If the latter
  is discrete (i.e. the nodes are multinomial random variables)
  each node has class \code{bn.fit.dnode} and contains the
  following elements:

  \itemize{

    \item \code{node}: the label of the node.
    \item \code{parents}: the labels of the parents of the node.
    \item \code{children}: the labels of the children of the node.
    \item \code{prob}: the conditional probability table of the node
      given its parents.

  }

  If on the other hand the network is continuous (i.e. the nodes
  are Gaussian random variables) each node has class \code{bn.fit.gnode}
  and contains the following elements:

  \itemize{

    \item \code{node}: the label of the node.
    \item \code{parents}: the labels of the parents of the node.
    \item \code{children}: the labels of the children of the node.
    \item \code{coefficients}: the linear regression coefficients of the
      parents against the node.
    \item \code{residuals}: the residuals of the linear regression,
      that is response minus fitted values.
    \item \code{fitted.values}: the fitted mean values of the linear
      regression.
    \item \code{sd}: the standard deviation of the residuals.

  }

  Furthermore, Bayesian network classifiers store the label of the training
  node in an additional attribute named \code{training}.
  
}
\author{ Marco Scutari }
\keyword{ classes }
