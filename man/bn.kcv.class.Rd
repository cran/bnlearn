\name{bn.kcv class}
\alias{bn.kcv class}
\alias{bn.kcv-class}
\alias{bn.kcv.list class}
\alias{bn.kcv.list-class}
\title{The bn.kcv class structure}
\description{

  The structure of an object of S3 class \code{bn.kcv} or \code{bn.kcv.list}.

}
\details{

  An object of class \code{bn.kcv.list} is a list whose elements are objects
  of class \code{bn.kcv}.

  An object of class \code{bn.kcv} is a list whose elements correspond to the
  iterations of a k-fold cross-validation. Each element contains the following
  objects:

  \itemize{

    \item \code{test}: an integer vector, the indexes of the observations
      used as a test set.
    \item \code{fitted}: an object of class \code{bn.fit}, the Bayesian network
      fitted from the training set.
    \item \code{loss}: the value of the loss function.

  }

  If the loss function requires to predict values from the test sets, each
  element also contains:

  \itemize{

    \item \code{predicted}: a factor or a numeric vector, the predicted values
      for the target node in the test set.
    \item \code{observed}: a factor or a numeric vector, the observed values
      for the target node in the test set.

  }

  In addition, an object of class \code{bn.kcv} has the following attributes:

  \itemize{

    \item \code{loss}: a character string, the label of the loss function.
    \item \code{mean}: the mean of the values of the loss function computed in
      the \code{k} iterations of the cross-validation.
    \item \code{bn}: either a character string (the label of the learning
      algorithm to be applied to the training data in each iteration) or an
      object of class \code{bn} (a fixed network structure).

  }

}
\author{ Marco Scutari }
\keyword{ classes }
