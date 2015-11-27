\name{cpquery}
\alias{cpquery}
\alias{cpdist}
\alias{mutilated}
\title{Perform conditional probability queries}
\description{

  Perform conditional probability queries (CPQs).

}
\usage{
cpquery(fitted, event, evidence, cluster = NULL, method = "ls", ...,
  debug = FALSE)
cpdist(fitted, nodes, evidence, cluster = NULL, method = "ls", ...,
  debug = FALSE)

mutilated(x, evidence)
}
\arguments{
  \item{fitted}{an object of class \code{bn.fit}.}
  \item{x}{an object of class \code{bn} or \code{bn.fit}.}
  \item{event, evidence}{see below.}
  \item{nodes}{a vector of character strings, the labels of the nodes whose
    conditional distribution we are interested in.}
  \item{cluster}{an optional cluster object from package \pkg{parallel}. See
    \code{\link{parallel integration}} for details and a simple example.}
  \item{method}{a character string, the method used to perform the conditional
    probability query. Currently only \emph{logic sampling} (\code{ls}, the
    default) and \emph{likelihood weighting} (\code{lw}) are implemented.}
  \item{\dots}{additional tuning parameters.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
}
\details{

  \code{cpquery} estimates the conditional probability of \code{event} given
  \code{evidence} using the method specified in the \code{method} argument.

  \code{cpdist} generates random observations conditional on the
  \code{evidence} using the method specified in the \code{method} argument.

  \code{mutilated} constructs the mutilated network used for sampling in
  \emph{likelihood weighting}.

  Note that both \code{cpquery} and \code{cpdist} are based on Monte Carlo
  particle filters, and therefore they may return slightly different values
  on different runs.

}
\section{Logic Sampling}{

  The \code{event} and \code{evidence} arguments must be two expressions
  describing the event of interest and the conditioning evidence in a format
  such that, if we denote with \code{data} the data set the network was learned
  from, \code{data[evidence, ]} and \code{data[event, ]} return the correct
  observations. If either \code{event} or \code{evidence} is set to \code{TRUE}
  an unconditional probability query is performed with respect to that argument.

  Three tuning parameters are available:

  \itemize{

    \item \code{n}: a positive integer number, the number of random observations
      to generate from \code{fitted}. The default value is 
      \code{5000 * log10(nparams.fitted(fitted))} for discrete and coditional
      Gaussian networks and \code{500 * nparams.fitted(fitted)} for Gaussian
      networks.

    \item \code{batch}: a positive integer number, the size of each batch of
     random observations. Defaults to \code{10^4}.

    \item \code{query.nodes}: a a vector of character strings, the labels of
      the nodes involved in \code{event} and \code{evidence}. Simple queries do
      not require to generate observations from all the nodes in the network,
      so \code{cpquery} and \code{cpdist} try to identify which nodes are used
      in \code{event} and \code{evidence} and reduce the network to their upper
      closure. \code{query.nodes} may be used to manually specify these nodes
      when automatic identification fails; there is no reason to use it 
      otherwise.

  }

  Note that the number of observations returned by \code{cpdist} is always
  smaller than \code{n}, because logic sampling is a form of rejection sampling.
  Therefore, only the obervations matching \code{evidence} (out of the \code{n}
  that are generated) are returned, and their number depends on the probability
  of \code{evidence}.

}
\section{Likelihood Weighting}{

  The \code{event} argument must be an expression describing the event of
  interest, as in logic sampling. The \code{evidence} argument must be a named
  list:
  \itemize{

    \item Each element corresponds to one node in the network and must contain
      the value that node will be set to when sampling.
    \item In the case of a continuous node, two values can also be provided. In
      that case, the value for that node will be sampled from a uniform
      distribution on the interval delimited by the specified values.
    \item In the case of a discrete or ordinal node, two or more values can also
      be provided. In that case, the value for that node will be sampled with
      uniform probability from the set of specified values.

  }

  If either \code{event} or \code{evidence} is set to \code{TRUE} an
  unconditional probability query is performed with respect to that argument.

  Tuning parameters are the same as for logic sampling: \code{n}, \code{batch}
  and \code{query.nodes}.

  Note that the observations returned by \code{cpdist} are generated from the
  mutilated network, and need to be weighted appropriately when computing
  summary statistics (for more details, see the references below).
  \code{cpquery} does that automatically when computing the final conditional
  probability. Also note that the \code{batch} argument is ignored in cpdist
  for speed and memory efficiency.

}
\value{

  \code{cpquery} returns a numeric value, the conditional probability of
  \code{event} conditional on \code{evidence}.

  \code{cpdist} returns a data frame containing the observations generated from
  the conditional distribution of the \code{nodes} conditional on
  \code{evidence}. The data frame has class \code{c("bn.cpdist", "data.frame")},
  and a \code{method} attribute storing the value of the \code{method} argument.
  In the case of likelihood weighting, the weights are also attached as an
  attribute called \code{weights}.

  \code{mutilated} returns a \code{bn} or \code{bn.fit} object, depending on the
  class of \code{x}.

}
\references{

  Koller D, Friedman N (2009). \emph{Probabilistic Graphical Models: Principles
      and Techniques}. MIT Press.

  Korb K, Nicholson AE (2010). \emph{Bayesian Artificial Intelligence}. Chapman
      & Hall/CRC, 2nd edition.

}
\examples{
## discrete Bayesian network (it is the same with ordinal nodes).
data(learning.test)
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
# do the same with likelihood weighting
cpquery(fitted, event = eval(parse(text = str2)),
  evidence = as.list(learning.test[2, -3]), method = "lw")
# conditional distribution of A given C == "c".
table(cpdist(fitted, "A", (C == "c")))

## Gaussian Bayesian network.
data(gaussian.test)
fitted = bn.fit(hc(gaussian.test), gaussian.test)
# the result should be around 0.04.
cpquery(fitted,
  event = ((A >= 0) & (A <= 1)) & ((B >= 0) & (B <= 3)),
  evidence = (C + D < 10))
}
\author{ Marco Scutari }
\keyword{multivariate}
