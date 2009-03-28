\name{bn-strength class}
\alias{bn-strength}
\title{

  The bn-strength class structure

}
\description{

  The structure of an object of the \code{bn-strength} S3 class.

}
\details{

  An object of class \code{bn-strength} is a data frame with some
  additional attributes:

  \itemize{

    \item \code{mode}: a character string, the criterion used to
      compute the strength coefficient. It can be equal to either
      \code{test} or \code{score}.
    \item \code{threshold}: a numeric value, the threshold used
      to determine if an arc is significant.

  }

  No method is defined for this class; therefore it can be
  manipulated as a standard data frame.

}
\author{ Marco Scutari }
\keyword{ internal }
