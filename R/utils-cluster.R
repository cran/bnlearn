# check the cluster is running.
isClusterRunning = function(cl) {

  tryCatch(any(unlist(clusterEvalQ(cl, TRUE))),
    error = function(err) { FALSE })

}#ISCLUSTERRUNNING

# check the status of the snow cluaster.
check.cluster = function(cluster) {

  if (is.null(cluster))
    return(TRUE)

  if (!(any(class(cluster) %in% supported.clusters)))
    stop("cluster is not a valid cluster object.")
  if (!(require(snow)))
    stop("Can't find required packages: snow")
  if (!isClusterRunning(cluster))
    stop("the cluster is stopped.")

}#CHECK.CLUSTER

# get the number of slaves.
nSlaves = function(cluster) {

  length(cluster)

}#NSLAVES
