\name{lizards}
\docType{data}
\alias{lizards}
\title{Lizards' perching behaviour data set}
\description{

  Real-world data set about the perching behaviour of two species of lizards in
  the South Bimini island, from Shoener (1968).

}
\usage{
data(lizards)
}
\format{

  The \code{lizards} data set contains the following variables:
  \itemize{

    \item \code{Species} (\emph{the species of the lizard}): a two-level factor
      with levels \code{Sagrei} and \code{Distichus}.
    \item \code{Height} (\emph{perch height}): a two-level factor with levels
      \code{high} (greater than 4.75 feet) and \code{low} (lesser or equal to
      4.75 feet).
    \item \code{Diameter} (\emph{perch diameter}): a two-level factor with
      levels \code{narrow} (greater than 4 inches) and \code{wide} (lesser or
      equal to 4 inches).

  }

}
\source{

  Edwards DI (2000). \emph{Introduction to Graphical Modelling}. Springer, 2nd
    edition.

  Fienberg SE (1980). \emph{The Analysis of Cross-Classified Categorical Data}.
    Springer, 2nd edition.

  Schoener TW (1968). "The Anolis Lizards of Bimini: Resource Partitioning in a
    Complex Fauna". \emph{Ecology}, \strong{49}(4):704--726.

}
\examples{
# load the data.
data(lizards)
# create and plot the network structure.
dag = model2network("[Species][Diameter|Species][Height|Species]")
\dontrun{graphviz.plot(dag, shape = "ellipse")}

# This data set is useful as it offers nominal values for
# the conditional mutual information and X^2 tests.
ci.test("Height", "Diameter", "Species", test = "mi", data = lizards)
ci.test("Height", "Diameter", "Species", test = "x2", data = lizards)
}
\keyword{datasets}
