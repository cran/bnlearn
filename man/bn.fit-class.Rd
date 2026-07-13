\name{bn.fit class}
\alias{bn.fit class}
\alias{bn.fit-class}
\alias{bn.fit.dnode}
\alias{bn.fit.gnode}
\title{The bn.fit class structure}
\description{

  The structure of an object of S3 class \code{bn.fit}.

}
\details{

  An object of class \code{bn.fit} is a list whose elements correspond to the
  nodes of the Bayesian network. If the network is discrete (the nodes are
  multinomial random variables), the object also has class \code{bn.fit.dnet}.
  Each node has class \code{bn.fit.dnode} and contains the following elements:

  \itemize{

    \item \code{node}: a character string, the label of the node.
    \item \code{parents}: a vector of character strings, the labels of the
      parents of the node.
    \item \code{children}: a vector of character strings, the labels of the
      children of the node.
    \item \code{prob}: a (multi)dimensional numeric table, the conditional
      probability table of the node given its parents.

  }

  Nodes encoding ordinal variables (ordered factors) have class
  \code{bn.fit.onode} and contain the same elements as \code{bn.fit.dnode}
  nodes. Networks with only ordinal nodes also have class \code{bn.fit.onet},
  while those containing both ordinal and multinomial nodes also have class
  \code{bn.fit.donet}.

  If, on the other hand, the network is continuous (that is, the nodes are
  Gaussian random variables), the object also has class \code{bn.fit.gnet}; each
  node has class \code{bn.fit.gnode} and contains the following elements:

  \itemize{

    \item \code{node}: a character string, the label of the node.
    \item \code{parents}: a vector of character strings, the labels of the
      parents of the node.
    \item \code{children}: a vector of character strings, the labels of the
      children of the node.
    \item \code{coefficients}: a numeric vector, the linear regression
      coefficients of the parents against the node.
    \item \code{residuals}: a numeric vector, the residuals of the linear
      regression.
    \item \code{fitted.values}: a numeric vector, the fitted mean values of
      the linear regression.
    \item \code{sd}: a numeric value, the standard deviation of the residuals
      (that is, the standard error).

  }

  Hybrid (conditional linear Gaussian) networks also have class
  \code{bn.fit.gnet}. Gaussian nodes have class \code{bn.fit.gnode}, discrete
  nodes have class \code{bn.fit.dnode} and conditional Gaussian nodes have class
  \code{bn.fit.cgnode}. Each node contains the following elements:

  \itemize{

    \item \code{node}: a character string, the label of the node.
    \item \code{parents}: a vector of character strings, the labels of the
      parents of the node.
    \item \code{children}: a vector of character strings, the labels of the
      children of the node.
    \item \code{dparents}: an integer vector, the indexes of the discrete
      parents in \code{parents}.
    \item \code{gparents}: an integer vector, the indexes of the continuous
      parents in \code{parents}.
    \item \code{dlevels}: a list containing the levels of the discrete
      parents in \code{parents}.
    \item \code{coefficients}: a numeric matrix, the linear regression
      coefficients of the continuous parents. Each column corresponds to
      a configuration of the discrete parents.
    \item \code{residuals}: a numeric vector, the residuals of the linear
      regression.
    \item \code{fitted.values}: a numeric vector, the fitted mean values of the
      linear regression.
    \item \code{configs}: an integer vector, the indexes of the configurations
      of the discrete parents.
    \item \code{sd}: a numeric vector, the standard deviation of the residuals
      (the standard error) for each configuration of the discrete parents.

  }

  Zero-Inflated (non-integer count) networks also have class \code{bn.fit.zihp}
  (hyper-Poisson) or \code{bn.fit.zinb} (negative binomial). Zero-inflated
  hyper-Poisson nodes have class \code{bn.fit.zihpnode} and contain the
  following elements:

  \itemize{

    \item \code{node}: a character string, the label of the node.
    \item \code{parents}: a vector of character strings, the labels of the
      parents of the node.
    \item \code{children}: a vector of character strings, the labels of the
      children of the node.
    \item \code{inflation}: a numeric vector, the coefficients of the
      zero-inflation (logistic) regression model associated with the parents.
    \item \code{intensity}: a numeric vector, the coefficients of the
      (log-linear) regression model for the Poisson intensity parameter
      associated with the parents.
    \item \code{dispersion}: a numeric value, the dispersion parameter that
      augments the Poisson (pseudo)counts into a hyper-Poisson distribution.

  }

  Zero-inflated negative binomial nodes have class \code{bn.fit.zinbnode} and
  contain the following elements:

  \itemize{

    \item \code{node}: a character string, the label of the node.
    \item \code{parents}: a vector of character strings, the labels of the
      parents of the node.
    \item \code{children}: a vector of character strings, the labels of the
      children of the node.
    \item \code{inflation}: a numeric vector, the coefficients of the
      zero-inflation (logistic) regression model associated with the parents.
    \item \code{prsucc}: a numeric vector, the coefficients of the
      (logistic) regression model for the probability of success parameter
      associated with the parents.
    \item \code{failures}: a numeric value, the number of failures.

  }

  In addition, Bayesian network classifiers store the label of the training node
  in an attribute named \code{training}.

  An object of class \code{bn.fit} may also have auxiliary class labels,
  depending on how it has been created. In particular:

  \itemize{

    \item \code{bn.naive}, identifying naive Bayes classifiers from
      \code{naive.bayes()};
    \item \code{bn.tan}, identifying tree-augmented naive Bayes (TAN)
      classifiers from \code{tree.bayes()};
    \item \code{bn.twin}, identifying twin networks from \code{twin()};
    \item \code{bn.ctf}, identifying counterfactual networks from \code{twin()}.

  }

}
\author{Marco Scutari}
\keyword{classes}
