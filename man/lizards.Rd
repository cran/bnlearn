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

  Edwards DI (2000). \emph{Introduction to Graphical Modelling}.
      Springer, 2nd edition.

  Fienberg SE (1980). \emph{The Analysis of Cross-Classified
      Categorical Data}. Springer, 2nd edition.

  Schoener TW (1968). "The Anolis Lizards of Bimini: Resource
      Partitioning in a Complex Fauna". \emph{Ecology},
      \strong{49}(4), 704-726.

}
\examples{

# load the data and build the correct network from the model string.
data(lizards)
res = empty.graph(names(lizards))
modelstring(res) = "[Species][Diameter|Species][Height|Species]"
plot(res)

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
# 	Mutual Information (discrete)
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
