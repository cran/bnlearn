\name{choose.direction}
\alias{choose.direction}
\title{ Try to infer the direction of an undirected arc }
\description{

  Check both possible directed arcs for existence,
  and choose the one with the lowest p-value.

}
\usage{
  choose.direction(x, arc, data, debug = FALSE)
}
\arguments{
  \item{x}{an object of class "bn".}
  \item{arc}{ a character string vector of length 2, the labels of
       two nodes of the graph.}
  \item{data}{a data frame, containing the data the Bayesian network
       was learned from.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output 
       is printed; otherwise the function is completely silent.}
}
\value{

  \code{choose.direction} returns invisibly an updated copy of \code{x}.

}
\examples{

data(learning.test)

res = gs(learning.test)
# the arc between E and F has no direction
plot(res)
res = choose.direction(res, c("E", "F"), learning.test, debug = TRUE)
# * testing E - F for direction.
#   > testing E -> F with conditioning set '  '.
#     > p-value is 2.174867e-197 .
#   > testing F -> E with conditioning set ' B '.
#     > p-value is 3.935648e-245 .
#   @ removing E -> F .
# * (re)building cached information about network structure.
#   > detecting neighbourhood, parents and children of node A .
#   > detecting neighbourhood, parents and children of node B .
#   > detecting neighbourhood, parents and children of node C .
#   > detecting neighbourhood, parents and children of node D .
#   > detecting neighbourhood, parents and children of node E .
#   > detecting neighbourhood, parents and children of node F .
# * detecting markov blanket for node A .
#   > neighbourhood is ' B D '.
#   > for child D getting ' A C '.
#   > raw markov blanket then is ' B D A C '.
#   > clean markov blanket then is ' B D C '.
# * detecting markov blanket for node B .
#   > neighbourhood is ' A E '.
#   > for child E getting ' B F '.
#   > raw markov blanket then is ' A E B F '.
#   > clean markov blanket then is ' A E F '.
# * detecting markov blanket for node C .
#   > neighbourhood is ' D '.
#   > for child D getting ' A C '.
#   > raw markov blanket then is ' D A C '.
#   > clean markov blanket then is ' D A '.
# * detecting markov blanket for node D .
#   > neighbourhood is ' A C '.
#   > raw markov blanket then is ' A C '.
#   > clean markov blanket then is ' A C '.
# * detecting markov blanket for node E .
#   > neighbourhood is ' B F '.
#   > raw markov blanket then is ' B F '.
#   > clean markov blanket then is ' B F '.
# * detecting markov blanket for node F .
#   > neighbourhood is ' E '.
#   > for child E getting ' B F '.
#   > raw markov blanket then is ' E B F '.
#   > clean markov blanket then is ' E B '.
plot(res, highlight = c("E", "F"))

}
\author{ Marco Scutari }
\keyword{htest}
