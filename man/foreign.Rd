\name{foreign files utilities}
\alias{read.bif}
\alias{write.bif}
\alias{read.dsc}
\alias{write.dsc}
\alias{read.net}
\alias{write.net}
\alias{write.dot}
\title{Read and write BIF, NET, DSC and DOT files}
\description{

  Read networks saved from other programs into \code{bn.fit} objects, and dump
  \code{bn} and \code{bn.fit} objects into files for other programs to read.

}
\usage{
# Old (non-XML) Bayesian Interchange format.
read.bif(file, debug = FALSE)
write.bif(file, fitted)

# Microsoft Interchange format.
read.dsc(file, debug = FALSE)
write.dsc(file, fitted)

# HUGIN flat network format.
read.net(file, debug = FALSE)
write.net(file, fitted)

# Graphviz DOT format.
write.dot(file, graph)
}
\arguments{
  \item{file}{a connection object or a character string.}
  \item{fitted}{an object of class \code{bn.fit}.}
  \item{graph}{an object of class \code{bn} or \code{bn.fit}.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
}
\note{

  All the networks present in the Bayesian Network Repository have associated
  BIF, DSC and NET files that can be imported with \code{read.bif}, 
  \code{read.dsc} and \code{read.net}.

  HUGIN can import and export NET files; Netica can read (but not write) DSC
  files; and Genie can read and write both DSC and NET files.

  DOT files can be read by Graphviz, Gephi and a variety of other programs.

  Please note that these functions work on a "best effort" basis, as the parsing
  of these formats have been implemented by reverse engineering the file format
  from publicly available examples.

}
\value{

  \code{read.bif}, \code{read.dsc} and \code{read.net} return an object of class
  \code{bn.fit}.

  \code{write.bif}, \code{write.dsc}, \code{write.net} and \code{write.dot}
  return \code{NULL} invisibly.

}
\references{

  Bayesian Network Repository, \url{http://www.bnlearn.com/bnrepository}.

  Genie, \url{http://genie.sis.pitt.edu}.

  HUGIN Expert, \url{http://www.hugin.com}.

  Netica, \url{http://www.norsys.com/netica.html}.

  Graphviz, \url{https://github.com/ellson/graphviz}.

  Gephi, \url{http://gephi.org}.

}
\author{ Marco Scutari }
\keyword{IO}
\keyword{file}
\keyword{utilities}
