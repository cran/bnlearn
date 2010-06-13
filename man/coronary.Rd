\name{coronary}
\docType{data}
\alias{coronary}
\title{Coronary Heart Disease data set}
\description{

  Probable risk factors for coronary trombosis, comprising
  data from 1841 men.

}
\usage{
data(coronary)
}
\format{

  The \code{coronary} data set contains the following 6 variables:
  \itemize{

    \item \code{Smoking} (\emph{smoking}): a two-level factor with
      levels \code{no} and \code{yes}.
    \item \code{M. Work} (\emph{strenuous mental work}): a two-level
      factor with levels \code{no} and \code{yes}.
    \item \code{P. Work} (\emph{strenuous physical work}): a two-level
      factor with levels \code{no} and \code{yes}.
    \item \code{Pressure} (\emph{systolic blood pressure}): a two-level
      factor with levels \code{<140} and \code{>140}.
    \item \code{Proteins} (\emph{ratio of beta and alpha lipoproteins}):
      a two-level factor with levels \code{<3} and \code{>3}.
    \item \code{Family} (\emph{family anamnesis of coronary heart
      disease}): a two-level factor with levels \code{neg} and \code{pos}.

  }

}
\source{

  Edwards DI (2000). \emph{Introduction to Graphical Modelling}.
      Springer, 2nd edition.

  Reinis Z, Pokorny J, Basika V, Tiserova J, Gorican K, Horakova D,
      Stuchlikova E, Havranek T, Hrabovsky F (1981). "Prognostic Significance
      of the Risk Profile in the Prevention of Coronary Heart Disease".
      Bratisl Lek Listy, \strong{76}, 137-150. Published on Bratislava
      Medical Journal, in Czech.

  Whittaker J (1990). \emph{Graphical Models in Applied Multivariate
      Statistics}. Wiley.

}
\examples{

# This is the undirected graphical model from Whittaker (1990).
data(coronary)
ug = empty.graph(names(coronary))
arcs(ug, ignore.cycles = TRUE) = matrix(
  c("Family", "M. Work", "M. Work", "Family",
    "M. Work", "P. Work", "P. Work", "M. Work",
    "M. Work", "Proteins", "Proteins", "M. Work",
    "M. Work", "Smoking", "Smoking", "M. Work",
    "P. Work", "Smoking", "Smoking", "P. Work",
    "P. Work", "Proteins", "Proteins", "P. Work",
    "Smoking", "Proteins", "Proteins", "Smoking",
    "Smoking", "Pressure", "Pressure", "Smoking",
    "Pressure", "Proteins", "Proteins", "Pressure"),
  ncol = 2, byrow = TRUE,
  dimnames = list(c(), c("from", "to")))

}
\keyword{datasets}
