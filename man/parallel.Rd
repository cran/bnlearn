\name{parallel integration}
\alias{parallel integration}
\title{bnlearn - snow/parallel package integration}
\description{

  How to use the \pkg{bnlearn} package with the parallel computing environment
  provided by the \pkg{snow} or \pkg{parallel} packages.

}
\section{Parallel computing for constraint-based algorithms}{
\preformatted{
# load parallel and bnlearn and rsprng.
> library(parallel)
> library(bnlearn)
> cl = makeCluster(2)
# check it works.
> clusterEvalQ(cl, runif(10))
[[1]]
 [1] 0.9245585 0.1876445 0.3371175 0.2267916 0.0392876 0.9085125
 [7] 0.9041345 0.7408525 0.1537343 0.9503611

[[2]]
 [1] 0.1932651 0.8218854 0.6087155 0.9037118 0.5257906 0.8737284
 [7] 0.5225114 0.8149691 0.1671706 0.6883363
# load the data.
> data(learning.test)
# call a learning function passing the cluster object (the
# return value of the previous makeCluster() call) as a
# parameter.
> res = gs(learning.test, cluster = cl)
# note that the number of test is evenly divided between
# the two nodes of the cluster.
> clusterEvalQ(cl, test.counter())
[[1]]
[1] 47

[[2]]
[1] 42
# a few tests are still executed by this process.
> (test.counter())
[1] 4
# stop the cluster.
> stopCluster(cl)
}
}
\author{ Marco Scutari }
\keyword{ documentation }
