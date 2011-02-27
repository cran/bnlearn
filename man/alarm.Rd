\name{alarm}
\docType{data}
\alias{alarm}
\title{ALARM Monitoring System (synthetic) data set}
\description{

  The ALARM ("A Logical Alarm Reduction Mechanism") is a Bayesian
  network designed to provide an alarm message system for patient
  monitoring.

}
\usage{
data(alarm)
}
\format{

  The \code{alarm} data set contains the following 37 variables:
  \itemize{

    \item \code{CVP} (\emph{central venous pressure}): a three-level
      factor with levels \code{LOW}, \code{NORMAL} and \code{HIGH}.
    \item \code{PCWP} (\emph{pulmonary capillary wedge pressure}):
      a three-level factor with levels \code{LOW}, \code{NORMAL} and
      \code{HIGH}.
    \item \code{HIST} (\emph{history}): a two-level factor with
      levels \code{TRUE} and \code{FALSE}.
    \item \code{TPR} (\emph{total peripheral resistance}): a three-level
      factor with levels \code{LOW}, \code{NORMAL} and \code{HIGH}.
    \item \code{BP} (\emph{blood pressure}): a three-level
      factor with levels \code{LOW}, \code{NORMAL} and \code{HIGH}.
    \item \code{CO} (\emph{cardiac output}): a three-level
      factor with levels \code{LOW}, \code{NORMAL} and \code{HIGH}.
    \item \code{HRBP} (\emph{heart rate / blood pressure}): a three-level
      factor with levels \code{LOW}, \code{NORMAL} and \code{HIGH}.
    \item \code{HREK} (\emph{heart rate measured by an EKG monitor}):
      a three-level factor with levels \code{LOW}, \code{NORMAL} and
      \code{HIGH}.
    \item \code{HRSA} (\emph{heart rate / oxygen saturation}): a three-level
      factor with levels \code{LOW}, \code{NORMAL} and \code{HIGH}.
    \item \code{PAP} (\emph{pulmonary artery pressure}): a three-level
      factor with levels \code{LOW}, \code{NORMAL} and \code{HIGH}.
    \item \code{SAO2} (\emph{arterial oxygen saturation}): a three-level
      factor with levels \code{LOW}, \code{NORMAL} and \code{HIGH}.
    \item \code{FIO2} (\emph{fraction of inspired oxygen}): a two-level
      factor with levels \code{LOW} and \code{NORMAL}.
    \item \code{PRSS} (\emph{breathing pressure}): a four-level factor
      with levels \code{ZERO}, \code{LOW}, \code{NORMAL} and \code{HIGH}.
    \item \code{ECO2} (\emph{expelled CO2}): a four-level factor with
      levels \code{ZERO}, \code{LOW}, \code{NORMAL} and \code{HIGH}.
    \item \code{MINV} (\emph{minimum volume}): a four-level factor with levels
      \code{ZERO}, \code{LOW}, \code{NORMAL} and \code{HIGH}.
    \item \code{MVS} (\emph{minimum volume set}): a three-level factor
      with levels \code{LOW}, \code{NORMAL} and \code{HIGH}.
    \item \code{HYP} (\emph{hypovolemia}): a two-level factor with
      levels \code{TRUE} and \code{FALSE}.
    \item \code{LVF} (\emph{left ventricular failure}): a two-level
      factor with levels \code{TRUE} and \code{FALSE}.
    \item \code{APL} (\emph{anaphylaxis}): a two-level factor with
      levels \code{TRUE} and \code{FALSE}.
    \item \code{ANES} (\emph{insufficient anesthesia/analgesia}): a
      two-level factor with levels \code{TRUE} and \code{FALSE}.
    \item \code{PMB} (\emph{pulmonary embolus}): a two-level factor
      with levels \code{TRUE} and \code{FALSE}.
    \item \code{INT} (\emph{intubation}): a three-level factor with
      levels \code{NORMAL}, \code{ESOPHAGEAL} and \code{ONESIDED}.
    \item \code{KINK} (\emph{kinked tube}): a two-level factor with
      levels \code{TRUE} and \code{FALSE}.
    \item \code{DISC} (\emph{disconnection}): a two-level factor with
      levels \code{TRUE} and \code{FALSE}.
    \item \code{LVV} (\emph{left ventricular end-diastolic volume}): a
      three-level factor with levels \code{LOW}, \code{NORMAL} and
      \code{HIGH}.
    \item \code{STKV} (\emph{stroke volume}): a three-level
      factor with levels \code{LOW}, \code{NORMAL} and \code{HIGH}.
    \item \code{CCHL} (\emph{catecholamine}): a two-level factor with
      levels \code{NORMAL} and \code{HIGH}.
    \item \code{ERLO} (\emph{error low output}): a two-level factor with
      levels \code{TRUE} and \code{FALSE}.
    \item \code{HR} (\emph{heart rate}): a three-level
      factor with levels \code{LOW}, \code{NORMAL} and \code{HIGH}.
    \item \code{ERCA} (\emph{electrocauter}): a two-level factor with
      levels \code{TRUE} and \code{FALSE}.
    \item \code{SHNT} (\emph{shunt}): a two-level factor with
      levels \code{NORMAL} and \code{HIGH}.
    \item \code{PVS} (\emph{pulmonary venous oxygen saturation}): a
      three-level factor with levels \code{LOW}, \code{NORMAL} and
      \code{HIGH}.
    \item \code{ACO2} (\emph{arterial CO2}): a three-level
      factor with levels \code{LOW}, \code{NORMAL} and \code{HIGH}.
    \item \code{VALV} (\emph{pulmonary alveoli ventilation}): a four-level
      factor with levels \code{ZERO}, \code{LOW}, \code{NORMAL} and
      \code{HIGH}.
    \item \code{VLNG} (\emph{lung ventilation}): a four-level factor with
      levels \code{ZERO}, \code{LOW}, \code{NORMAL} and \code{HIGH}.
    \item \code{VTUB} (\emph{ventilation tube}): a four-level factor with
      levels \code{ZERO}, \code{LOW}, \code{NORMAL} and \code{HIGH}.
    \item \code{VMCH} (\emph{ventilation machine}): a four-level factor with
      levels \code{ZERO}, \code{LOW}, \code{NORMAL} and \code{HIGH}.

  }

}
\source{

  Beinlich I, Suermondt HJ, Chavez RM, Cooper GF (1989). "The ALARM Monitoring
      System: A Case Study with Two Probabilistic Inference Techniques for
      Belief Networks." In "Proceedings of the 2nd European Conference on
      Artificial Intelligence in Medicine", pp. 247-256. Springer-Verlag.

  Elidan G (2001). "Bayesian Network Repository". \cr
      \url{http://www.cs.huji.ac.il/site/labs/compbio/Repository}.

}
\note{

  The R script to generate data from this network is shipped in the
  \file{network.scripts} directory of this package.

}
\examples{
# load the data and build the correct network from the model string.
data(alarm)
res = empty.graph(names(alarm))
modelstring(res) = paste("[HIST|LVF][CVP|LVV][PCWP|LVV][HYP][LVV|HYP:LVF]",
  "[LVF][STKV|HYP:LVF][ERLO][HRBP|ERLO:HR][HREK|ERCA:HR][ERCA]",
  "[HRSA|ERCA:HR][ANES][APL][TPR|APL][ECO2|ACO2:VLNG][KINK]",
  "[MINV|INT:VLNG][FIO2][PVS|FIO2:VALV][SAO2|PVS:SHNT][PAP|PMB][PMB]",
  "[SHNT|INT:PMB][INT][PRSS|INT:KINK:VTUB][DISC][MVS][VMCH|MVS]",
  "[VTUB|DISC:VMCH][VLNG|INT:KINK:VTUB][VALV|INT:VLNG][ACO2|VALV]",
  "[CCHL|ACO2:ANES:SAO2:TPR][HR|CCHL][CO|HR:STKV][BP|CO:TPR]", sep = "")
\dontrun{
# there are too many nodes for plot(), use graphviz.plot().
graphviz.plot(res)}}
\keyword{datasets}
