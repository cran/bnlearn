\name{foreign files utilities}
\alias{read.bif}
\alias{write.bif}
\alias{read.dsc}
\alias{write.dsc}
\alias{read.net}
\alias{write.net}
\title{ Read and write BIF, NET and DSC files }
\description{

  Read networks saved from other programs into \code{bn.fit} objects, and
  dump \code{bn.fit} objects into files for other programs to read.

}
\usage{
# Old (non-XML) Bayesian Interchange files.
read.bif(file, debug = FALSE)
write.bif(file, fitted)

# Microsoft Interchange files.
read.dsc(file, debug = FALSE)
write.dsc(file, fitted)

# HUGIN flat network format.
read.net(file, debug = FALSE)
write.net(file, fitted)
}
\arguments{
  \item{file}{a connection object or a character string.}
  \item{fitted}{an object of class \code{bn.fit}.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output
    is printed; otherwise the function is completely silent.}
}
\note{

  Most of the networks present in the Bayesian Network Repository
  have an associated BIF file that can be imported with \code{read.bif}.

  HUGIN can import and export NET files; Netica can read (but not write)
  DSC files; and Genie can read and write both DSC and NET files.

  Please note that these functions work on a "best effort" basis, as
  the parsing of these formats have been implementing by reverse
  engineering the file format from publicly available examples.

}
\value{

  \code{read.bif}, \code{read.dsc} and \code{read.net} return an
    object of class \code{bn.fit}.

  \code{write.bif}, \code{write.dsc} and \code{write.net} return
    \code{NULL} invisibly.

}
\references{

  Elidan G (2001). "Bayesian Network Repository". \cr
      \url{http://www.cs.huji.ac.il/site/labs/compbio/Repository}.

  Genie, \url{http://genie.sis.pitt.edu}.

  HUGIN Expert, \url{http://www.hugin.com}.

  Netica, \url{http://www.norsys.com/netica.html}.

}
\author{ Marco Scutari }
\keyword{IO}
\keyword{file}
\keyword{utilities}
