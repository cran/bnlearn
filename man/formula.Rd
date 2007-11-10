\name{model string tools}
\alias{modelstring}
\alias{model2network}
\alias{as.character.bn}
\alias{as.bn}
\alias{as.bn.character}
\title{ Build a model string from a Bayesian network and vice versa }
\description{

  Build a model string from a Bayesian network and vice versa.

}
\usage{

  modelstring(x)
  model2network(string, debug = FALSE)

  \method{as.character}{bn}(x, ...)
  \method{as.bn}{character}(x, debug = FALSE)

}
\arguments{
  \item{x}{an object of class "bn".}
  \item{string}{a character string describing the Bayesian network.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output 
       is printed; otherwise the function is completely silent.}
  \item{...}{extra arguments from the generic method (currently ignored).}
}
\details{

  The strings returned by \code{modelstring} have the same format as
  the ones returned by the \code{modelstring} function in package
  \pkg{deal}; network structures may be easily exported to and imported
  from that package (via the \code{model2network} function). 

}
\value{

  \code{model2network} and \code{as.bn.character} return an object of 
  class \code{bn}; \code{modelstring} and \code{as.character.bn} return
  a character string.

}
\examples{

data(learning.test)
res = set.arc(gs(learning.test), "A", "B")

modelstring(res)
# [1] "[F][C][A][B|A][D|A:C][E|B:F]"
res2 = model2network(modelstring(res))
compare(res, res2)
# [1] TRUE

}
\author{ Marco Scutari }
\keyword{utilities}
