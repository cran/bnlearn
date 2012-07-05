\name{cpquery}
\alias{cpquery}
\alias{cpdist}
\title{ Perform conditional probability queries }
\description{

  Perform conditional probability queries (CPQs).

}
\usage{
cpquery(fitted, event, evidence, cluster = NULL, method = "ls", ...,
  debug = FALSE)
cpdist(fitted, nodes, evidence, cluster = NULL, method = "ls", ...,
  debug = FALSE)
}
\arguments{
  \item{fitted}{an object of class \code{bn.fit}.}
  \item{event, evidence}{see below.}
  \item{nodes}{a vector of character strings, the labels of the nodes whose
    conditional distribution we are interested in.}
  \item{cluster}{an optional cluster object from package \pkg{snow}.
      See \code{\link{snow integration}} for details and a simple
      example.}
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

   Note that the number of observations returned by \code{cpdist} is always 
   smaller than \code{n}, because logic sampling is a form of rejection 
   sampling. Therefore, conly the obervations matching \code{evidence} (out
   of the \code{n} that are generated) are returned, and their number depends
   on the probability of \code{evidence}.

}
\value{

  \code{cpquery} returns a numeric value, the conditional probability of \code{event}
  conditional on \code{evidence}.

  \code{cpudist} returns a data frame containing the observations generated from
  the conditional distribution of the \code{nodes} conditional on \code{evidence}.

}
\references{

  Koller D, Friedman N (2009). \emph{Probabilistic Graphical Models:
      Principles and Techniques}. MIT Press.

  Korb K, Nicholson AE (2010). \emph{Bayesian Artificial Intelligence}.
      Chapman & Hall/CRC, 2nd edition.

}
\examples{
fitted = bn.fit(hc(learning.test), learning.test)
# the result should be around 0.025.
cpquery(fitted, (B == "b"), (A == "a"))
# for a single observation, predict the value of a single
# variable conditional on the others.
var = names(learning.test)
obs = 2
str = paste("(", names(learning.test)[-3], "=='",
        sapply(learning.test[obs,-3], as.character), "')", 
        sep = "", collapse = " & ")
str
str2 = paste("(", names(learning.test)[3], "=='",
         as.character(learning.test[obs, 3]), "')", sep = "")
str2
cpquery(fitted, eval(parse(text = str2)), eval(parse(text = str)))
# conditional distribution of A given C == "c". 
table(cpdist(fitted, "A", (C == "c")))
}
\author{ Marco Scutari }
\keyword{multivariate}
