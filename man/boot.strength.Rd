\name{boot.strength}
\alias{boot.strength}
\title{ Bootstrap arc strength and direction }
\description{

  Use nonparametric bootstrap to assess arc strength
  and direction.

}
\usage{
boot.strength(data, R = 200, m = nrow(data),
  algorithm, algorithm.args = list(), debug = FALSE)
}
\arguments{
  \item{data}{a data frame containing the variables in the model.}
  \item{R}{a positive integer, the number of bootstrap replicates.}
  \item{m}{a positive integer, the size of each bootstrap replicate.}
  \item{algorithm}{a character string, the learning algorithm to be
    applied to the bootstrap replicates. Possible values are \code{gs},
    \code{iamb}, \code{fast.iamb}, \code{inter.iamb}, \code{mmpc}
    and \code{hc}. See \code{\link{bnlearn-package}} and the
    documentation of each algorithm for details.}
  \item{algorithm.args}{a list of extra arguments to be passed to
    the learning algorithm.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output
       is printed; otherwise the function is completely silent.}
}
\value{

  A 4-column data frame (very similar to an object of class
  \code{bn.strength}) with an entry for each possible arc in
  the network and the following columns:

    \item{from, to}{the nodes incident on the arc.}
    \item{strength}{the strength of the arc, computed as the probability
      of observing an arc between \code{from} and \code{to} in the bootstrap
      replicates, regardless of its direction.}
    \item{direction}{the confidence in the direction of the arc, computed
      as the probability of that particular direction in the bootstrap
      replicates conditional on the presence of an arc between \code{from}
      and \code{to} (again regardless of its direction).}

}
\references{

  Imoto S, Kim SY, Shimodaira H, Aburatani S, Tashiro K, Kuhara S,
      Miyano S (2002). "Bootstrap Analysis of Gene Networks Based on
      Bayesian Networks and Nonparametric Regression". \emph{Genome
      Informatics}, \strong{13}, 369-370.

}
\author{ Marco Scutari }
\seealso{\code{\link{arc.strength}}, \code{\link{bn.boot}}. }
\keyword{nonparametric}
\keyword{multivariate}
