\name{graphviz.chart}
\alias{graphviz.chart}
\title{Plotting networks with probability bars}
\description{

  Plot a Bayesian network as a graph whose nodes are barplots representing
  the marginal probability distributions of the corresponding variables.
  Requires the \pkg{Rgraphviz} and \pkg{gRain} packages.

}
\usage{
graphviz.chart(x, type = "barchart", layout = "dot", draw.levels = TRUE,
  grid = FALSE, scale = c(0.75, 1.1), col = "black", bg = "transparent",
  text.col = "black", bar.col = "black", strip.bg = bg, main = NULL,
  sub = NULL)
}
\arguments{
  \item{x}{an object of class \code{bn.fit} representing a discrete Bayesian
    network.}
  \item{type}{a character string, the type of graph used to plot the probability
    distributions in the nodes. Possible values are \code{barchart}, \code{dotplot}
    and \code{barprob} (a barchart with probabilities printed over the bars).}
  \item{layout}{a character string, the layout argument that will be passed to
    \pkg{Rgraphviz}. Possible values are \code{dots}, \code{neato},
    \code{twopi}, \code{circo} and \code{fdp}. See \pkg{Rgraphviz} documentation
    for details.}
  \item{draw.levels}{a boolean value, whether to print the labels of the levels
    of each variable.}
  \item{grid}{a boolean value, whether to draw to a reference grid for the
    probability distributions. If \code{grid} is \code{TRUE}, a vertical grid is
    drawn at probabilities \code{c(0, 0.25,} \code{0.50, 0.75)}. If \code{grid}
    is a numeric vector, a verical grid is drawn at the specified probabilities.}
  \item{scale}{a vector of two positive numbers, used by \pkg{Rgraphviz} to
    determine the size and the aspect ratio of the nodes.}
  \item{col, bg, text.col, bar.col, strip.bg}{the colours of the node border,
    of the barchart background, of the text, of the bars and of the strip
    background.}
  \item{main}{a character string, the main title of the graph. It's plotted at
    the top of the graph.}
  \item{sub}{a character string, a subtitle which is plotted at the bottom of
    the graph.}
}
\value{

  \code{graphviz.chart()} invisibly returns \code{NULL}.

}
\examples{
\dontrun{
modelstring = paste("[HIST|LVF][CVP|LVV][PCWP|LVV][HYP][LVV|HYP:LVF][LVF]",
  "[STKV|HYP:LVF][ERLO][HRBP|ERLO:HR][HREK|ERCA:HR][ERCA][HRSA|ERCA:HR][ANES]",
  "[APL][TPR|APL][ECO2|ACO2:VLNG][KINK][MINV|INT:VLNG][FIO2][PVS|FIO2:VALV]",
  "[SAO2|PVS:SHNT][PAP|PMB][PMB][SHNT|INT:PMB][INT][PRSS|INT:KINK:VTUB][DISC]",
  "[MVS][VMCH|MVS][VTUB|DISC:VMCH][VLNG|INT:KINK:VTUB][VALV|INT:VLNG][ACO2|VALV]",
  "[CCHL|ACO2:ANES:SAO2:TPR][HR|CCHL][CO|HR:STKV][BP|CO:TPR]", sep = "")
dag = model2network(modelstring)
fitted = bn.fit(dag, alarm)

# Netica style.
graphviz.chart(fitted, grid = TRUE, bg = "beige", bar.col = "black")
# Hugin style.
graphviz.chart(fitted, type = "barprob", grid = TRUE, bar.col = "green",
  strip.bg = "lightyellow")
# GeNIe style.
graphviz.chart(fitted, col = "darkblue", bg = "azure", bar.col = "darkblue")
# personal favourites.
graphviz.chart(fitted, type = "barprob", grid = TRUE, bar.col = "darkgreen",
  strip.bg = "lightskyblue")
graphviz.chart(fitted, type = "barprob", grid = TRUE, bar.col = "gold",
  strip.bg = "lightskyblue")
# dot-plot version.
graphviz.chart(fitted, type = "dotplot")
}
}
\author{Marco Scutari}
\keyword{plots}
\keyword{interfaces to other packages}
