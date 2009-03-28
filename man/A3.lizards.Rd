\name{lizards}
\docType{data}
\alias{lizards}
\title{Lizards' perching behaviour data set}
\description{

  Real-world data set about the perching behaviour of two
  species of lizards in the South Bimini island, from
  Shoener (1968).

}
\usage{
  data(lizards)
}
\format{

  The \code{lizards} data set contains the following variables:
  \itemize{

    \item \code{Species} (\emph{the species of the lizard}): a
      two-level factor with levels \code{Sagrei} and
      \code{Distichus}.
    \item \code{Height} (\emph{perch height}): a two-level
      factor with levels \code{high} (greater than 4.75 feet)
      and \code{low} (lesser or equal to 4.75 feet).
    \item \code{Diameter} (\emph{perch diameter}): a
      two-level factor with levels \code{narrow} (greater
      than 4 inches) and \code{wide} (lesser or equal to 4
      inches).

  }

}
\source{

  D. Edwards. Introduction to Graphical Modelling.
    Springer, 2000.

  S. E. Fienberg. The Analysis of Cross-Classified Categorical
    Data. MIT Press, 1980.

  T. W. Schoener (1968). The anolis lizards of Bimini: Resource
    partitioning in a complex fauna. Ecology, 49: 704 -- 726.

}
\examples{

## The graphical model for this data set is:
#       Height - Species - Diameter

data(lizards)
table(lizards[, c(3,2,1)])
# , , Species = Sagrei
#
#       Diameter
# Height narrow wide
#   high     86   35
#   low      32   11
#
# , , Species = Distichus
#
#       Diameter
# Height narrow wide
#   high     73   70
#   low      61   41

\dontrun{
# This data set is useful as it offers nominal values for
# the conditional mutual information and X^2 tests.

attach(lizards)
ci.test(Height, Diameter, Species, test = "mi")
#
# 	mutual information (discrete)
#
# data:  Height ~ Diameter | Species
# mi = 2.0256, df = 2, p-value = 0.3632
# alternative hypothesis: true value is greater than 0
ci.test(Height, Diameter, Species, test = "x2")
#
# 	Pearson's X^2
#
# data:  Height ~ Diameter | Species
# x2 = 2.0174, df = 2, p-value = 0.3647
# alternative hypothesis: true value is greater than 0
}

}
\keyword{datasets}
