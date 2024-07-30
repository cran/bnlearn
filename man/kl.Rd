\name{information theoretic quantities}
\alias{H}
\alias{KL}
\title{Compute the distance between two fitted Bayesian networks}
\description{

  Compute Shannon's entropy of a fitted Bayesian network and the
  Kullback-Leibler divergence between two fitted Bayesian networks.

}
\usage{
H(P)
KL(P, Q)
}
\arguments{
  \item{P, Q}{objects of class \code{bn.fit}.}
}
\note{

  Note that in the case of Gaussian and conditional Gaussian netwoks the
  divergence can be negative. Regardless of the type of network, if at least one
  of the two networks is singular the divergence can be infinite.

  If any of the parameters of the two networks are \code{NA}s, the divergence
  will also be \code{NA}.

}
\value{

  \code{H()} and \code{KL()} return a single numeric value.

}
\examples{
\dontrun{
# discrete networks
dag = model2network("[A][C][F][B|A][D|A:C][E|B:F]")
fitted1 = bn.fit(dag, learning.test, method = "mle")
fitted2 = bn.fit(dag, learning.test, method = "bayes", iss = 20)

H(fitted1)
H(fitted2)

KL(fitted1, fitted1)
KL(fitted2, fitted2)
KL(fitted1, fitted2)
}

# continuous, singular networks.
dag = model2network("[A][B][E][G][C|A:B][D|B][F|A:D:E:G]")
singular = fitted1 = bn.fit(dag, gaussian.test)
singular$A = list(coef = coef(fitted1[["A"]]) + runif(1), sd = 0)

H(singular)
H(fitted1)

KL(singular, fitted1)
KL(fitted1, singular)
}
\author{Marco Scutari}
\keyword{parameter learning}
\keyword{inference}
