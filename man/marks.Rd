\name{marks}
\docType{data}
\alias{marks}
\title{Examination marks data set}
\description{

  Examination marks of 88 students on five different
  topics, from Mardia (1979).

}
\usage{
data(marks)
}
\format{

  The \code{marks} data set contains the following variables,
  one for each topic in the examination:
  \itemize{

    \item \code{MECH} (\emph{mechanics})
    \item \code{VECT} (\emph{vectors})
    \item \code{ALG} (\emph{algebra})
    \item \code{ANL} (\emph{analysis})
    \item \code{STAT} (\emph{statistics})

  }

  All are measured on the same scale (0-100).

}
\source{

  Edwards DI (2000). \emph{Introduction to Graphical Modelling}.
      Springer, 2nd edition.

  Mardia KV, Kent JT, Bibby JM (1979). \emph{Multivariate Analysis}.
      Academic Press.

  Whittaker J (1990). \emph{Graphical Models in Applied Multivariate
      Statistics}. Wiley.


}
\examples{

# This is the undirected graphical model from Edwards (2000).
data(marks)
ug = empty.graph(names(marks))
arcs(ug, ignore.cycles = TRUE) = matrix(
  c("MECH", "VECT", "MECH", "ALG", "VECT", "MECH", "VECT", "ALG",
    "ALG",  "MECH", "ALG", "VECT", "ALG",  "ANL", "ALG",  "STAT",
    "ANL",  "ALG", "ANL",  "STAT", "STAT", "ALG", "STAT", "ANL"),
  ncol = 2, byrow = TRUE,
  dimnames = list(c(), c("from", "to")))

}
\keyword{datasets}
