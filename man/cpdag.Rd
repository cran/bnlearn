\name{cpadg}
\alias{cpdag}
\title{ Find the equivalence class of a Bayesian nework  }
\description{

  Find the equivalence class of a Bayesian nework.

}
\usage{
  cpdag(x, debug = FALSE)
}
\arguments{
  \item{x}{an object of class \code{bn}.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output
       is printed; otherwise the function is completely silent.}
}
\value{

  An object of class \code{bn}, representing the equivalence class.
  See \code{\link{bn-class}} for details.

}
\references{

  Pearl J (2000). \emph{Causality: Models, Reasning and Inference}. 
      Cambridge University Press.

}
\examples{
data(learning.test)
cpdag(gs(learning.test))
}
\author{ Marco Scutari }
\keyword{utilities}
\keyword{multivariate}
