# check whether the cluster is running.
isClusterRunning = function(cl) {

  tryCatch(any(unlist(parallel::clusterEvalQ(cl, TRUE))),
    error = function(err) { FALSE })

}#ISCLUSTERRUNNING

# check the status of the snow/parallel cluster.
check.cluster = function(cluster) {

  if (is.null(cluster))
    return(TRUE)
  if (!is(cluster, supported.clusters))
    stop("cluster is not a valid cluster object.")
  if (!requireNamespace("parallel"))
    stop("this function requires the parallel package.")
  if (!isClusterRunning(cluster))
    stop("the cluster is stopped.")

}#CHECK.CLUSTER

# get the number of slaves.
nSlaves = function(cluster) {

  length(cluster)

}#NSLAVES

slaves.setup = function(cluster) {

  # set the test counter in all the cluster nodes.
  parallel::clusterEvalQ(cluster, library(bnlearn))
  parallel::clusterEvalQ(cluster, reset.test.counter())

}#SLAVE.SETUP

# smart parLapply() that falls back to standard lapply().
smartLapply = function(cl, ...) {

  if (is.null(cl))
    lapply(...)
  else
    parallel::parLapply(cl = cl, ...)

}#SMARTLAPPLY

