\name{compare}
\alias{compare}
\alias{all.equal.bn}
\alias{shd}
\alias{hamming}
\alias{graphviz.compare}
\title{Compare two or more different Bayesian networks}
\description{

  Compare two different Bayesian networks; compute their Structural Hamming
  Distance (SHD) or the Hamming distance between their skeletons. Or
  graphically compare them by plotting them side by side,

}
\usage{
compare(target, current, arcs = FALSE)
\method{all.equal}{bn}(target, current, ...)

shd(learned, true, wlbl = FALSE, debug = FALSE)
hamming(learned, true, debug = FALSE)

graphviz.compare(x, ..., groups, layout = "dot", shape = "circle", main = NULL,
  sub = NULL, diff = "from-first", diff.args = list())
}
\arguments{
  \item{target, learned}{an object of class \code{bn}.}
  \item{current, true}{another object of class \code{bn}.}
  \item{\dots}{extra arguments from the generic method (for \code{all.equal()},
    currently ignored); or a set of one or more objects of class \code{bn}
    (for \code{graphviz.compare}).}
  \item{wlbl}{a boolean value. If \code{TRUE} arcs whose directions have been
    fixed by a whitelist or a by blacklist are preserved when constructing
    the CPDAGs of \code{learned} and \code{true}.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
  \item{arcs}{a boolean value. See below.}
  \item{x}{an object of class \code{bn}.}
  \item{groups}{a list of character vectors, representing groups of node labels
    of nodes that should be plotted close to each other.}
  \item{layout}{a character string, the layout argument that will be passed to
    \pkg{Rgraphviz}. Possible values are \code{dots}, \code{neato},
    \code{twopi}, \code{circo} and \code{fdp}. See \pkg{Rgraphviz}
    documentation for details.}
  \item{shape}{a character string, the shape of the nodes. Can be \code{circle},
    \code{ellipse} or \code{rectangle}.}
  \item{main}{a vector of character strings, one for each network. They
    are plotted at the top of the corresponding figure(s).}
  \item{sub}{a vector of character strings, the subtitles that are plotted at
    the bottom of the corresponding figure(s).}
  \item{diff}{a character string, the label of the method used to compare and
    format the figure(s) created by \code{graphviz.compare()}. The default value
    is \code{from-first}, se below for details.}
  \item{diff.args}{a list of optional arguments to control the formatting of
    the figure(s) created by \code{graphviz.compare()}. See below for details.}
}
\details{

  \code{graphviz.compare()} can visualize differences between graphs in various
  way depending on the value of the \code{diff} and \code{diff.args} arguments:
  \itemize{

    \item \code{none}: differences are not highlighted.
    \item \code{from-first}: the first \code{bn} object, \code{x}, is taken as
      the reference network. All the other networks, passed via the \code{\dots}
      argument, are compared to that first network and their true positive,
      false positive, false negative arcs relative to that first network are
      highlighted. Colours, line types and line widths for each category of
      arcs can be specified as the elements of a list via the \code{diff.args}
      argument, with names \code{tp.col}, \code{tp.lty}, \code{tp.lwd},
      \code{fp.col}, \code{fp.lty}, \code{fp.lwd}, \code{fn.col}, \code{fn.lty},
      \code{tp.lwd}. In addition, it is possible not to plot the reference
      network at all by setting \code{show.first} to \code{FALSE}.

  }

  Regardless of the visualization, the nodes are arranged to be in the same
  position for all the networks to make it easier to compare them.

}
\note{

  Note that SHD, as defined in the reference, is defined on CPDAGs; therefore
  \code{cpdag()} is called on both \code{learned} and \code{true} before computing
  the distance.

}
\value{

  \code{compare()} returns a list containing the number of true positives
  (\code{tp}, the number of arcs in \code{current} also present in
  \code{target}), of false positives (\code{fp}, the number of arcs in
  \code{current} not present in \code{target}) and of false negatives
  (\code{fn}, the number of arcs not in \code{current} but present in
  \code{target}) if \code{arcs} is \code{FALSE}; or the corresponding arc sets
  if \code{arcs} is \code{TRUE}.

  \code{all.equal()} returns either \code{TRUE} or a character string describing
  the differences between \code{target} and \code{current}.

  \code{shd()} and \code{hamming()} return a non-negative integer number.

  \code{graphviz.compare()} plots one or more figures and returns invisibly a
  list containing the \code{graph} objects generated from the networks that were
  passed as arguments (in the same order). They can be further modified using
  the \pkg{graph} and \pkg{Rgraphviz} packages.
}
\examples{
data(learning.test)

e1 = model2network("[A][B][C|A:B][D|B][E|C][F|A:E]")
e2 = model2network("[A][B][C|A:B][D|B][E|C:F][F|A]")
shd(e2, e1, debug = TRUE)
unlist(compare(e1,e2))
compare(target = e1, current = e2, arcs = TRUE)
\dontrun{graphviz.compare(e1, e2, diff = "none")}
}
\references{

  Tsamardinos I, Brown LE, Aliferis CF (2006). "The Max-Min Hill-Climbing
    Bayesian Network Structure Learning Algorithm". \emph{Machine Learning},
    \strong{65}(1):31--78.

}
\author{Marco Scutari}
\keyword{graphs}
\keyword{plots}
\keyword{interfaces to other packages}
