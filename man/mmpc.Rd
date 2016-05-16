\name{local discovery algorithms}
\alias{local discovery algorithms}
\alias{aracne}
\alias{chow.liu}
\title{Local discovery structure learning algorithms}
\description{

  ARACNE and Chow-Liu learn simple graphs structures from data using pairwise
  mutual information coefficients.

}
\usage{
aracne(x, whitelist = NULL, blacklist = NULL, mi = NULL, debug = FALSE)
chow.liu(x, whitelist = NULL, blacklist = NULL, mi = NULL, debug = FALSE)
}
\arguments{
  \item{x}{a data frame containing the variables in the model.}
  \item{whitelist}{a data frame with two columns (optionally labeled "from"
    and "to"), containing a set of arcs to be included in the graph.}
  \item{blacklist}{a data frame with two columns (optionally labeled "from"
    and "to"), containing a set of arcs not to be included in the graph.}
  \item{mi}{a character string, the estimator used for the pairwise (i.e.
    unconditional) mutual information coefficients in the ARACNE and Chow-Liu
    algorithms. Possible values are \code{mi} (discrete mutual information)
    and \code{mi-g} (Gaussian mutual information).}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
}
\value{

  An object of class \code{bn}. See \code{\link{bn-class}} for details.

}
\references{

  Margolin AA, Nemenman I, Basso K, Wiggins C, Stolovitzky G, Dalla Favera R,
      Califano A (2006). "ARACNE: An Algorithm for the Reconstruction of
      Gene Regulatory Networks in a Mammalian Cellular Context". \emph{BMC
      Bioinformatics}, \strong{7}(Suppl 1):S7.

}
\author{ Marco Scutari }
\seealso{\link{constraint-based algorithms}, \link{score-based algorithms},
  \link{hybrid algorithms}.}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}
