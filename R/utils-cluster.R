# check the cluster is running.
isClusterRunning = function(cl) {

  tryCatch(any(clusterEvalQ(cl, TRUE)),
    error = function(err) { FALSE })

}#ISCLUSTERRUNNING


