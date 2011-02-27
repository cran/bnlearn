\name{hailfinder}
\docType{data}
\alias{hailfinder}
\title{The HailFinder weather forecast system (synthetic) data set}
\description{

  Hailfinder is a Bayesian network designed to forecast severe
  summer hail in northeastern Colorado.

}
\usage{
data(hailfinder)
}
\format{

  The \code{hailfinder} data set contains the following 56 variables:
  \itemize{

    \item \code{N07muVerMo} (\emph{10.7mu vertical motion}): a four-level factor with levels \code{StrongUp}, \code{WeakUp}, \code{Neutral} and \code{Down}.
    \item \code{SubjVertMo} (\emph{subjective judgment of vertical motion}): a four-level factor with levels \code{StrongUp}, \code{WeakUp}, \code{Neutral} and \code{Down}.
    \item \code{QGVertMotion} (\emph{quasigeostrophic vertical motion}): a four-level factor with levels \code{StrongUp}, \code{WeakUp}, \code{Neutral} and \code{Down}.
    \item \code{CombVerMo} (\emph{combined vertical motion}): a four-level factor with levels \code{StrongUp}, \code{WeakUp}, \code{Neutral} and \code{Down}.
    \item \code{AreaMesoALS} (\emph{area of meso-alpha}): a four-level factor with levels \code{StrongUp}, \code{WeakUp}, \code{Neutral} and \code{Down}.
    \item \code{SatContMoist} (\emph{satellite contribution to moisture}): a four-level factor with levels \code{VeryWet}, \code{Wet}, \code{Neutral} and \code{Dry}.
    \item \code{RaoContMoist} (\emph{reading at the forecast center for moisture}): a four-level factor with levels \code{VeryWet}, \code{Wet}, \code{Neutral} and \code{Dry}.
    \item \code{CombMoisture} (\emph{combined moisture}): a four-level factor with levels \code{VeryWet}, \code{Wet}, \code{Neutral} and \code{Dry}.
    \item \code{AreaMoDryAir} (\emph{area of moisture and adry air}): a four-level factor with levels \code{VeryWet}, \code{Wet}, \code{Neutral} and \code{Dry}.
    \item \code{VISCloudCov} (\emph{visible cloud cover}): a three-level factor with levels \code{Cloudy}, \code{PC} and \code{Clear}.
    \item \code{IRCloudCover} (\emph{infrared cloud cover}): a three-level factor with levels \code{Cloudy}, \code{PC} and \code{Clear}.
    \item \code{CombClouds} (\emph{combined cloud cover}): a three-level factor with levels \code{Cloudy}, \code{PC} and \code{Clear}.
    \item \code{CldShadeOth} (\emph{cloud shading, other}): a three-level factor with levels \code{Cloudy}, \code{PC} and \code{Clear}.
    \item \code{AMInstabMt} (\emph{AM instability in the mountains}): a three-level factor with levels \code{None}, \code{Weak} and \code{Strong}.
    \item \code{InsInMt} (\emph{instability in the mountains}): a three-level factor with levels \code{None}, \code{Weak} and \code{Strong}.
    \item \code{WndHodograph} (\emph{wind hodograph}): a four-level factor with levels \code{DCVZFavor}, \code{StrongWest}, \code{Westerly} and \code{Other}.
    \item \code{OutflowFrMt} (\emph{outflow from mountains}): a three-level factor with levels \code{None}, \code{Weak} and \code{Strong}.
    \item \code{MorningBound} (\emph{morning boundaries}): a three-level factor with levels \code{None}, \code{Weak} and \code{Strong}.
    \item \code{Boundaries} (\emph{boundaries}): a three-level factor with levels \code{None}, \code{Weak} and \code{Strong}.
    \item \code{CldShadeConv} (\emph{cloud shading, convection}): a three-level factor with levels \code{None}, \code{Some} and \code{Marked}.
    \item \code{CompPlFcst} (\emph{composite plains forecast}): a three-level factor with levels \code{IncCapDecIns}, \code{LittleChange} and \code{DecCapIncIns}.
    \item \code{CapChange} (\emph{capping change}): a three-level factor with levels \code{Decreasing}, \code{LittleChange} and \code{Increasing}.
    \item \code{LoLevMoistAd} (\emph{low-level moisture advection}): a four-level factor with levels \code{StrongPos}, \code{WeakPos}, \code{Neutral} and \code{Negative}.
    \item \code{InsChange} (\emph{instability change}):  three-level factor with levels \code{Decreasing}, \code{LittleChange} and \code{Increasing}.
    \item \code{MountainFcst} (\emph{mountains (region 1) forecast}): a three-level factor with levels \code{XNIL}, \code{SIG} and \code{SVR}.
    \item \code{Date} (\emph{date}): a six-level factor with levels \code{May15_Jun14}, \code{Jun15_Jul1}, \code{Jul2_Jul15}, \code{Jul16_Aug10}, \code{Aug11_Aug20} and \code{Aug20_Sep15}.
    \item \code{Scenario} (\emph{scenario}): an eleven-level factor with levels \code{A}, \code{B}, \code{C}, \code{D}, \code{E}, \code{F}, \code{G}, \code{H}, \code{I}, \code{J} and \code{K}.
    \item \code{ScenRelAMCIN} (\emph{scenario relevant to AM convective inhibition}): a two-level factor with levels \code{AB} and \code{CThruK}.
    \item \code{MorningCIN} (\emph{morning convective inhibition}): a four-level factor with levels \code{None}, \code{PartInhibit}, \code{Stifling} and \code{TotalInhibit}.
    \item \code{AMCINInScen} (\emph{AM convective inhibition in scenario}): a three-level factor with levels \code{LessThanAve}, \code{Average} and \code{MoreThanAve}.
    \item \code{CapInScen} (\emph{capping withing scenario}): a three-level factor with levels \code{LessThanAve}, \code{Average} and \code{MoreThanAve}.
    \item \code{ScenRelAMIns} (\emph{scenario relevant to AM instability}): a six-level factor with levels \code{ABI}, \code{CDEJ}, \code{F}, \code{G}, \code{H} and \code{K}.
    \item \code{LIfr12ZDENSd} (\emph{LI from 12Z DEN sounding}): a four-level factor with levels \code{LIGt0}, \code{N1GtLIGt_4}, \code{N5GtLIGt_8} and \code{LILt_8}.
    \item \code{AMDewptCalPl} (\emph{AM dewpoint calculations, plains}): a three-level factor with levels \code{Instability}, \code{Neutral} and \code{Stability}.
    \item \code{AMInsWliScen} (\emph{AM instability within scenario}): a three-level factor with levels \code{LessUnstable}, \code{Average} and \code{MoreUnstable}.
    \item \code{InsSclInScen} (\emph{instability scaling within scenario}): a three-level factor with levels \code{LessUnstable}, \code{Average} and \code{MoreUnstable}.
    \item \code{ScenRel34} (\emph{scenario relevant to regions 2/3/4}): a five-level factor with levels \code{ACEFK}, \code{B}, \code{D}, \code{GJ} and \code{HI}.
    \item \code{LatestCIN} (\emph{latest convective inhibition}): a four-level factor with levels \code{None}, \code{PartInhibit}, \code{Stifling} and \code{TotalInhibit}.
    \item \code{LLIW} (\emph{LLIW severe weather index}): a four-level factor with levels \code{Unfavorable}, \code{Weak}, \code{Moderate} and \code{Strong}.
    \item \code{CurPropConv} (\emph{current propensity to convection}): a four-level factor with levels \code{None}, \code{Slight}, \code{Moderate} and \code{Strong}.
    \item \code{ScnRelPlFcst} (\emph{scenario relevant to plains forecast}): an eleven-level factor with levels \code{A}, \code{B}, \code{C}, \code{D}, \code{E}, \code{F}, \code{G}, \code{H}, \code{I}, \code{J} and \code{K}.
    \item \code{PlainsFcst} (\emph{plains forecast}): a three-level factor with levels \code{XNIL}, \code{SIG} and \code{SVR}.
    \item \code{N34StarFcst} (\emph{regions 2/3/4 forecast}): a three-level factor with levels \code{XNIL}, \code{SIG} and \code{SVR}.
    \item \code{R5Fcst} (\emph{region 5 forecast}): a three-level factor with levels \code{XNIL}, \code{SIG} and \code{SVR}.
    \item \code{Dewpoints} (\emph{dewpoints}): a seven-level factor with levels \code{LowEverywhere}, \code{LowAtStation}, \code{LowSHighN}, \code{LowNHighS}, \code{LowMtsHighPl}, \code{HighEverywher}, \code{Other}.
    \item \code{LowLLapse} (\emph{low-level lapse rate}): a four-level factor with levels \code{CloseToDryAd}, \code{Steep}, \code{ModerateOrLe} and \code{Stable}.
    \item \code{MeanRH} (\emph{mean relative humidity}): a three-level factor with levels \code{VeryMoist}, \code{Average} and \code{Dry}.
    \item \code{MidLLapse} (\emph{mid-level lapse rate}): a three-level factor with levels \code{CloseToDryAd}, \code{Steep} and \code{ModerateOrLe}.
    \item \code{MvmtFeatures} (\emph{movement of features}): a four-level factor with levels \code{StrongFront}, \code{MarkedUpper}, \code{OtherRapid} and \code{NoMajor}.
    \item \code{RHRatio} (\emph{realtive humidity ratio}): a three-level factor with levels \code{MoistMDryL}, \code{DryMMoistL} and \code{other}.
    \item \code{SfcWndShfDis} (\emph{surface wind shifts and discontinuities}): a seven-level factor with levels \code{DenvCyclone}, \code{E_W_N}, \code{E_W_S}, \code{MovigFtorOt}, \code{DryLine}, \code{None} and \code{Other}.
    \item \code{SynForcng} (\emph{synoptic forcing}): a five-level factor with levels \code{SigNegative}, \code{NegToPos}, \code{SigPositive}, \code{PosToNeg} and \code{LittleChange}.
    \item \code{TempDis} (\emph{temperature discontinuities}): a four-level factor with levels \code{QStationary}, \code{Moving}, \code{None}, \code{Other}.
    \item \code{WindAloft} (\emph{wind aloft}): a four-level factor with levels \code{LV}, \code{SWQuad}, \code{NWQuad}, \code{AllElse}.
    \item \code{WindFieldMt} (\emph{wind fields, mountains}): a two-level factor with levels \code{Westerly} and \code{LVorOther}.
    \item \code{WindFieldPln} (\emph{wind fields, plains}): a six-level factor with levels \code{LV}, \code{DenvCyclone}, \code{LongAnticyc}, \code{E_NE}, \code{SEquad} and \code{WidespdDnsl}.

  }

}
\source{

  Abramson B, Brown J, Edwards W, Murphy A, Winkler RL (1996).
      "Hailfinder: A Bayesian system for forecasting severe weather".
      \emph{International Journal of Forecasting}, \strong{12}(1),
      57-71.

  Elidan G (2001). "Bayesian Network Repository". \cr
      \url{http://www.cs.huji.ac.il/site/labs/compbio/Repository}.

}
\note{

  The R script to generate data from this network is shipped in the
  \file{network.scripts} directory of this package.

}
\examples{
# load the data and build the correct network from the model string.
data(hailfinder)
res = empty.graph(names(hailfinder))
modelstring(res) = paste("[N07muVerMo][SubjVertMo][QGVertMotion]",
  "[SatContMoist][RaoContMoist][VISCloudCov][IRCloudCover][AMInstabMt]",
  "[WndHodograph][MorningBound][LoLevMoistAd][Date][MorningCIN]",
  "[LIfr12ZDENSd][AMDewptCalPl][LatestCIN][LLIW]",
  "[CombVerMo|N07muVerMo:SubjVertMo:QGVertMotion]",
  "[CombMoisture|SatContMoist:RaoContMoist]",
  "[CombClouds|VISCloudCov:IRCloudCover][Scenario|Date]",
  "[CurPropConv|LatestCIN:LLIW][AreaMesoALS|CombVerMo]",
  "[ScenRelAMCIN|Scenario][ScenRelAMIns|Scenario][ScenRel34|Scenario]",
  "[ScnRelPlFcst|Scenario][Dewpoints|Scenario][LowLLapse|Scenario]",
  "[MeanRH|Scenario][MidLLapse|Scenario][MvmtFeatures|Scenario]",
  "[RHRatio|Scenario][SfcWndShfDis|Scenario][SynForcng|Scenario]",
  "[TempDis|Scenario][WindAloft|Scenario][WindFieldMt|Scenario]",
  "[WindFieldPln|Scenario][AreaMoDryAir|AreaMesoALS:CombMoisture]",
  "[AMCINInScen|ScenRelAMCIN:MorningCIN]",
  "[AMInsWliScen|ScenRelAMIns:LIfr12ZDENSd:AMDewptCalPl]",
  "[CldShadeOth|AreaMesoALS:AreaMoDryAir:CombClouds]",
  "[InsInMt|CldShadeOth:AMInstabMt][OutflowFrMt|InsInMt:WndHodograph]",
  "[CldShadeConv|InsInMt:WndHodograph][MountainFcst|InsInMt]",
  "[Boundaries|WndHodograph:OutflowFrMt:MorningBound]",
  "[CompPlFcst|AreaMesoALS:CldShadeOth:Boundaries:CldShadeConv]",
  "[CapChange|CompPlFcst][InsChange|CompPlFcst:LoLevMoistAd]",
  "[CapInScen|CapChange:AMCINInScen]",
  "[InsSclInScen|InsChange:AMInsWliScen]",
  "[PlainsFcst|CapInScen:InsSclInScen:CurPropConv:ScnRelPlFcst]",
  "[N34StarFcst|ScenRel34:PlainsFcst][R5Fcst|MountainFcst:N34StarFcst]",
  sep = "")
\dontrun{
# there are too many nodes for plot(), use graphviz.plot().
graphviz.plot(res)}}
\keyword{datasets}
