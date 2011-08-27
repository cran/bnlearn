\name{bn.kcv class}
\alias{bn.kcv class}
\alias{bn.kcv-class}
\title{

  The bn.kcv class structure

}
\description{

  The structure of an object of the \code{bn.kcv} S3 class.

}
\details{

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

  If the loss function is some form of prediction error, each element also
  contains:

  \itemize{

    \item \code{predicted}: a factor, the predicted values for the test set.

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
