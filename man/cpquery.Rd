\name{cpquery}
\alias{cpquery}
\title{ Perform conditional probability queries }
\description{

  Perform conditional probability queries (CPQs).

}
\usage{
cpquery(fitted, event, evidence, method = "ls", ..., debug = FALSE)
}
\arguments{
  \item{fitted}{an object of class \code{bn.fit}.}
  \item{event, evidence}{see below.}
  \item{method}{a character string, the method used to perform the conditional
    probability query. Currently only \emph{Logic Sampling} is implemented.}
  \item{\dots}{additional tuning parameters.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output
       is printed; otherwise the function is completely silent.}
}
\section{Logic Sampling}{

   The \code{event} and \code{evidence} arguments must be two expressions
   describing the event of interest and the conditioning evidence in a
   format such that, if we denote with \code{data} the data set the network
   was learned from, \code{data[evidence, ]} and \code{data[event, ]}
   return the correct observations.
   If either parameter is equal to \code{TRUE} an unconditional probability
   query is performed.

   Two tuning parameters are available:

   \itemize{

     \item \code{n}: a positive integer number, the number of random observations
       to generate from \code{fitted}. Defaults to \code{5000 * nparams(fitted)}.

     \item \code{batch}: a positive integer number, the size of each batch
       of random observations. Defaults to \code{10^4}.

   }

}
\value{

  A numeric value, the conditional probability of event \code{conditional}
  on \code{evidence}.

}
\references{

  Korb K, Nicholson AE (2003). \emph{Bayesian Artificial Intelligence}.
      Chapman & Hall/CRC.

}
\examples{
fitted = bn.fit(hc(learning.test), learning.test)
cpquery(fitted, (B == "b"), (A == "a"))
# the result should be around 0.025.
}
\author{ Marco Scutari }
\keyword{multivariate}
