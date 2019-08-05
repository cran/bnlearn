\name{BF}
\alias{BF}
\title{Bayes factor between two network structures}
\description{

  Compute the Bayes factor between the structures of two Bayesian networks..

}
\usage{
BF(num, den, data, score, ..., log = TRUE)
}
\arguments{
  \item{num, den}{two objects of class \code{bn}, corresponding to the numerator
    and the denominator models in the Bayes factor.}
  \item{data}{a data frame containing the data to be used to compute the Bayes
    factor.}
  \item{score}{a character string, the label of a posterior network score.
    If none is specified, the default score is the \emph{Bayesian Dirichlet
    equivalent} score (\code{bde}) for discrete networks and the \emph{Bayesian
    Gaussian score} (\code{bge}) for Gaussian networks. Other kinds of Bayesian
    networks are not currently supported.}
  \item{\dots}{extra tuning arguments for the posterior scores. See
    \code{\link{score}} for details.}
  \item{log}{a boolean value. If \code{TRUE} the Bayes factor is given as
    log(BF).}
}
\note{

  The Bayes factor for two network structures, by definition, is the ratio of
  the respective marginal likelihoods which is equivalent to the ration of
  the corresponding posterior probabilities if we assume the \code{uniform}
  prior over all possible DAGs. However, note that it is possible to specify
  different priors using the \dQuote{\code{...}} arguments of \code{BF()}; in
  that case the value returned by the function will not be the classic Bayes
  factor.

}
\value{

  A single numeric value, the Bayes factor of the two network structures
  \code{num} and \code{den}.

}
\examples{
data(learning.test)

dag1 = model2network("[A][B][F][C|B][E|B][D|A:B:C]")
dag2 = model2network("[A][C][B|A][D|A][E|D][F|A:C:E]")
BF(dag1, dag2, learning.test, score = "bds", iss = 1)
}
\seealso{\code{\link{score}}, \code{\link{compare}},
  \code{\link{bf.strength}}.}
\author{Marco Scutari}
\keyword{network scores}
\keyword{structure learning}
