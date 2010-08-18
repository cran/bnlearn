\name{snow integration}
\alias{snow integration}
\title{

  bnlearn - snow package integration

}
\description{

  How to use the \pkg{bnlearn} package with the parallel computing
  environment provided by the \pkg{snow} package.

}
\section{Parallel computing for constraint-based algorithms}{
\preformatted{
# load snow, bnlearn and rsprng (for parallel random number
# generation, just in case it's needed); start LAM/MPI via
# lamboot if using an MPI cluster.
> library(snow)
> library(bnlearn)
> library(rsprng)
# initialize the cluster ("socket" and "PVM" clusters are fine, too).
> cl <- makeCluster(2, type = "MPI")
Loading required package: Rmpi
        2 slaves are spawned successfully. 0 failed.
> clusterSetupSPRNG(cl)
# load the data.
> data(learning.test)
# call a learning function passing the cluster object (the
# return value of the previous makeCluster() call) as a
# parameter.
> res = gs(learning.test, cluster = cl)
# note that the number of test is evenly divided between
# the two nodes of the cluster.
> clusterEvalQ(cl, .test.counter)
[[1]]
[1] 47

[[2]]
[1] 42
# a few tests are still executed by this process.
> .test.counter
[1] 4
# stop the cluster.
> stopCluster(cl)
[1] 1
}
}
\author{ Marco Scutari }
\keyword{ documentation }
