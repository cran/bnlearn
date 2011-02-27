\name{insurance}
\docType{data}
\alias{insurance}
\title{Insurance evaluation network (synthetic) data set}
\description{

  Insurance is a network for evaluating car insurance risks.

}
\usage{
data(insurance)
}
\format{

  The \code{insurance} data set contains the following 27 variables:
  \itemize{

    \item \code{GoodStudent} (\emph{good student}): a two-level factor with levels \code{False} and \code{True}.
    \item \code{Age} (\emph{age}): a three-level factor with levels \code{Adolescent}, \code{Adult} and \code{Senior}.
    \item \code{SocioEcon} (\emph{socio-economic status}): a four-level factor with levels \code{Prole}, \code{Middle}, \code{UpperMiddle} and \code{Wealthy}.
    \item \code{RiskAversion} (\emph{risk aversion}): a four-level factor with levels \code{Psychopath}, \code{Adventurous}, \code{Normal} and \code{Cautious}.
    \item \code{VehicleYear} (\emph{vehicle age}): a two-level factor with levels \code{Current} and \code{older}.
    \item \code{ThisCarDam} (\emph{damage to this car}): a four-level factor with levels \code{None}, \code{Mild}, \code{Moderate} and \code{Severe}.
    \item \code{RuggedAuto} (\emph{ruggedness of the car}): a three-level factor with levels \code{EggShell}, \code{Football} and \code{Tank}.
    \item \code{Accident} (\emph{severity of the accident}): a four-level factor with levels \code{None}, \code{Mild}, \code{Moderate} and \code{Severe}.
    \item \code{MakeModel} (\emph{car's model}): a five-level factor with levels \code{SportsCar}, \code{Economy}, \code{FamilySedan}, \code{Luxury} and \code{SuperLuxury}.
    \item \code{DrivQuality} (\emph{driving quality}): a three-level factor with levels \code{Poor}, \code{Normal} and \code{Excellent}.
    \item \code{Mileage} (\emph{mileage}): a four-level factor with levels \code{FiveThou}, \code{TwentyThou}, \code{FiftyThou} and \code{Domino}.
    \item \code{Antilock} (\emph{ABS}): a two-level factor with levels \code{False} and \code{True}.
    \item \code{DrivingSkill} (\emph{driving skill}): a three-level factor with levels \code{SubStandard}, \code{Normal} and \code{Expert}.
    \item \code{SeniorTrain} (\emph{senior training}): a two-level factor with levels \code{False} and \code{True}.
    \item \code{ThisCarCost} (\emph{costs for the insured car}): a four-level factor with levels \code{Thousand}, \code{TenThou}, \code{HundredThou} and \code{Million}.
    \item \code{Theft} (\emph{theft}): a two-level factor with levels \code{False} and \code{True}.
    \item \code{CarValue} (\emph{value of the car}): a five-level factor with levels \code{FiveThou}, \code{TenThou}, \code{TwentyThou}, \code{FiftyThou} and \code{Million}.
    \item \code{HomeBase} (\emph{neighbourhood type}): a four-level factor with levels \code{Secure}, \code{City}, \code{Suburb} and \code{Rural}.
    \item \code{AntiTheft} (\emph{anti-theft system}): a two-level factor with levels \code{False} and \code{True}.
    \item \code{PropCost} (\emph{ratio of the cost for the two cars}): a four-level factor with levels \code{Thousand}, \code{TenThou}, \code{HundredThou} and \code{Million}.
    \item \code{OtherCarCost} (\emph{costs for the other car}): a four-level factor with levels \code{Thousand}, \code{TenThou}, \code{HundredThou} and \code{Million}.
    \item \code{OtherCar} (\emph{other cars involved in the accident}): a two-level factor with levels \code{False} and \code{True}.
    \item \code{MedCost} (\emph{cost of the medical treatment}): a four-level factor with levels \code{Thousand}, \code{TenThou}, \code{HundredThou} and \code{Million}.
    \item \code{Cushioning} (\emph{cushioning}): a four-level factor with levels \code{Poor}, \code{Fair}, \code{Good} and \code{Excellent}.
    \item \code{Airbag} (\emph{airbag}): a two-level factor with levels \code{False} and \code{True}.
    \item \code{ILiCost} (\emph{inspection cost}): a four-level factor with levels \code{Thousand}, \code{TenThou}, \code{HundredThou} and \code{Million}.
    \item \code{DrivHist} (\emph{driving history}): a three-level factor with levels \code{Zero}, \code{One} and \code{Many}.

  }

}
\source{

  Binder J, Koller D, Russell S, Kanazawa K (1997). "Adaptive Probabilistic
      Networks with Hidden Variables". \emph{Machine Learning},
      \strong{29}(2-3), 213-244.

  Elidan G (2001). "Bayesian Network Repository". \cr
      \url{http://www.cs.huji.ac.il/site/labs/compbio/Repository}.

}
\note{

  The R script to generate data from this network is shipped in the
  \file{network.scripts} directory of this package.

}
\examples{
# load the data and build the correct network from the model string.
data(insurance)
res = empty.graph(names(insurance))
modelstring(res) = paste("[Age][Mileage][SocioEcon|Age]",
  "[GoodStudent|Age:SocioEcon][RiskAversion|Age:SocioEcon]",
  "[OtherCar|SocioEcon][VehicleYear|SocioEcon:RiskAversion]",
  "[MakeModel|SocioEcon:RiskAversion][SeniorTrain|Age:RiskAversion]",
  "[HomeBase|SocioEcon:RiskAversion][AntiTheft|SocioEcon:RiskAversion]",
  "[RuggedAuto|VehicleYear:MakeModel][Antilock|VehicleYear:MakeModel]",
  "[DrivingSkill|Age:SeniorTrain][CarValue|VehicleYear:MakeModel:Mileage]",
  "[Airbag|VehicleYear:MakeModel][DrivQuality|RiskAversion:DrivingSkill]",
  "[Theft|CarValue:HomeBase:AntiTheft][Cushioning|RuggedAuto:Airbag]",
  "[DrivHist|RiskAversion:DrivingSkill]",
  "[Accident|DrivQuality:Mileage:Antilock]",
  "[ThisCarDam|RuggedAuto:Accident][OtherCarCost|RuggedAuto:Accident]",
  "[MedCost|Age:Accident:Cushioning][ILiCost|Accident]",
  "[ThisCarCost|ThisCarDam:Theft:CarValue]",
  "[PropCost|ThisCarCost:OtherCarCost]", sep = "")
\dontrun{
# there are too many nodes for plot(), use graphviz.plot().
graphviz.plot(res)}}
\keyword{datasets}
