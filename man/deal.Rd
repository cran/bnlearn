\name{deal integration}
\alias{deal integration}
\title{

  bnlearn - deal package integration

}
\description{

  How to use the \pkg{bnlearn} package with the Bayesian
  network learning methods  provided by the \pkg{deal} package.

}
\section{Export a bn object to deal}{
\preformatted{
# load the bnlearn package.
> library(bnlearn)
> data(learning.test)
# learn the network structure.
> res = hc(learning.test)
> modelstring(res)
[1] "[A][C][F][B|A][D|A:C][E|B:F]"
# load the deal package.
> library(deal)

Attaching package: 'deal'


        The following object(s) are masked from package:bnlearn :

         modelstring,
         nodes,
         score

> bnlearn::node.ordering(res)
[1] "A" "C" "F" "B" "D" "E"
# create an empty network object.
> net = deal::network(learning.test[, node.ordering(res)])
# convert the bn object via its string representation.
> net = deal::as.network(bnlearn::modelstring(res), net)
# the network is the same, modulo some differences due to the
# partial ordering of the nodes.
> deal::modelstring(net)
[1] "[A][C][F][B|A][D|A:C][E|F:B]"
> bnlearn::modelstring(res)
[1] "[A][C][F][B|A][D|A:C][E|B:F]"
}
}
\section{Import a network structure from deal}{
\preformatted{
res2 = bnlearn::model2network(deal::modelstring(net))
}
}
\author{ Marco Scutari }
\keyword{ documentation }
