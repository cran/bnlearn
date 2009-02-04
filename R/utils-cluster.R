# check the cluster is running.
isClusterRunning = function(cl) {

  tryCatch(any(unlist(clusterEvalQ(cl, TRUE))),
    error = function(err) { FALSE })

}#ISCLUSTERRUNNING


