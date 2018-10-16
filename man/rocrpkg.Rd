\name{ROCR integration}
\alias{ROCR integration}
\alias{as.prediction}
\alias{as.prediction.bn.strength}
\title{Generating a prediction object for ROCR}
\description{

  Evaluate structure learning accuracy with \pkg{ROCR}. This function views the
  arcs in a \code{bn.strength} object as a set of predictions and the arcs in a
  \code{true} reference graph as a set of labels, and produces a \code{prediction}
  object from the \pkg{ROCR} package. This facilitates evaluation of structure
  learning with traditional machine learning metrics such as ROC curves and AUC.

}
\usage{
\method{as.prediction}{bn.strength}(x, true, ..., consider.direction = TRUE)
}
\arguments{
  \item{x}{an object of class \code{bn.strength} returned by \code{boot.strength()},
    representing learning results targeting the object of class \code{bn}
    specified by the \code{true} argument.}
  \item{true}{an object of class \code{bn}, the target of structure learning.}
  \item{\dots}{additional arguments, currently ignored.}
  \item{consider.direction}{a boolean value. If \code{TRUE} an arc's prediction
    value is set to the product of its \code{strength} and \code{direction}
    values in \code{x} (interpreted as the probability an arc is both present
    and has the specified direction). If \code{FALSE} the arc's prediction value
    is set to its \code{strength} value.}
}
\details{

  One way of evaluating the overall performance of a network structure learning
  algorithm is to evaluate how well it detects individual arcs.
  \code{as.prediction()} takes each pair of nodes in a ground truth network and
  labels them with a \code{1} if an arc exists between them and \code{0} if not.
  It uses the arc presence probabilities in a \code{bn.strength} object
  returned by \code{boot.strength()} as the predictions.

}
\value{

  An object of class \code{prediction} from the \pkg{ROCR} package.

}
\examples{
\dontrun{
library(ROCR)

modelstring = paste0("[HIST|LVF][CVP|LVV][PCWP|LVV][HYP][LVV|HYP:LVF][LVF]",
  "[STKV|HYP:LVF][ERLO][HRBP|ERLO:HR][HREK|ERCA:HR][ERCA][HRSA|ERCA:HR][ANES]",
  "[APL][TPR|APL][ECO2|ACO2:VLNG][KINK][MINV|INT:VLNG][FIO2][PVS|FIO2:VALV]",
  "[SAO2|PVS:SHNT][PAP|PMB][PMB][SHNT|INT:PMB][INT][PRSS|INT:KINK:VTUB][DISC]",
  "[MVS][VMCH|MVS][VTUB|DISC:VMCH][VLNG|INT:KINK:VTUB][VALV|INT:VLNG][ACO2|VALV]",
  "[CCHL|ACO2:ANES:SAO2:TPR][HR|CCHL][CO|HR:STKV][BP|CO:TPR]")
true.dag = model2network(modelstring)
strength = boot.strength(alarm, R = 200, m = 30, algorithm = "hc")
pred = as.prediction(strength, true.dag)
perf = performance(pred, "tpr", "fpr")
plot(perf, main = "Arc Detection")
performance(pred, "auc")
}
}
\author{Robert Ness}
\keyword{interfaces to other packages}
