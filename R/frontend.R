
gs = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL, 
    test = "mi", alpha = 0.05, debug = FALSE, optimized = TRUE, 
    strict = TRUE, direction = FALSE) {

  bnlearn(x = x, cluster = cluster, whitelist = whitelist, 
    blacklist = blacklist, test = test, alpha = alpha, debug = debug, 
    optimized = optimized, strict = strict, direction = direction)

}#GS

iamb = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL, 
    test = "mi", alpha = 0.05, debug = FALSE, optimized = TRUE, 
    strict = TRUE, direction = FALSE) {

  bnlearn(x = x, cluster = cluster, whitelist = whitelist, 
    blacklist = blacklist, test = test, alpha = alpha, method = "iamb", 
    debug = debug, optimized = optimized, strict = strict, 
    direction = direction)

}#IAMB

fast.iamb = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL, 
    test = "mi", alpha = 0.05, debug = FALSE, optimized = TRUE, 
    strict = TRUE, direction = FALSE) {
           
  bnlearn(x = x, cluster = cluster, whitelist = whitelist, 
    blacklist = blacklist, test = test, alpha = alpha, 
    method = "fast-iamb", debug = debug, optimized = optimized, 
    strict = strict, direction = direction)

}#FAST.IAMB

inter.iamb = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL, 
    test = "mi", alpha = 0.05, debug = FALSE, optimized = TRUE, 
    strict = TRUE, direction = FALSE) {
           
  bnlearn(x = x, cluster = cluster, whitelist = whitelist, 
    blacklist = blacklist, test = test, alpha = alpha, 
    method = "inter-iamb", debug = debug, optimized = optimized, 
    strict = strict, direction = direction)

}#INTER.IAMB

available.discrete.tests = c("mh", "mi", "fmi")
available.continuous.tests = c("cor", "zf")
available.tests = c(available.discrete.tests, available.continuous.tests)

bnlearn = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL, 
    test = "mh", alpha = 0.05, method = "gs", debug = FALSE, optimized = TRUE,
    strict = TRUE, direction = FALSE) {

  assign("test.counter", 0, envir = .GlobalEnv) 

  res = NULL
  available.methods = c("gs", "iamb", "fast-iamb", "inter-iamb")
  supported.clusters = c("MPIcluster", "PVMcluster","SOCKcluster")
  cluster.aware = FALSE

  # check the data.
  if(!is.data.frame(x)) 
    stop("x must be a data frame.")
  # check the test.
  if (!(test %in% available.tests))
    stop(paste("valid values for test are:", 
           paste(available.tests, collapse = " ")))
  # check the algorithm.
  if (!(method %in% available.methods))
    stop(paste("valid values for method are:", 
           paste(available.methods, collapse = " ")))
  # check if it's the right test for the data (discrete, continuous).
  data.class = sapply(x, class)
  if (any(data.class != "factor") && (test %in% available.discrete.tests))
    stop(paste("test '", test, "' may be used with discrete data only.", sep = ""))
  if (any(data.class != "numeric") && (test %in% available.continuous.tests))
    stop(paste("test '", test, "' may be used with continuous data only.", sep = ""))
  # check debug.
  if (!is.logical(debug) || is.na(debug))
    stop("debug must be a logical value (TRUE/FALSE).")
  # check strict.
  if (!is.logical(strict) || is.na(strict))
    stop("strict must be a logical value (TRUE/FALSE).")
  # check optimized.
  if (!is.logical(optimized) || is.na(optimized))
    stop("optimized must be a logical value (TRUE/FALSE).")
  # check direction.
  if (!is.logical(direction) || is.na(direction))
    stop("direction must be a logical value (TRUE/FALSE).")
  # check alpha.
  if (!is.numeric(alpha) || (alpha > 1) || (alpha < 0))
    stop("alpha must be a numerical value in [0,1].")

  # check cluster

  if (!is.null(cluster)) {

    if (!(any(class(cluster) %in% supported.clusters)))
      stop("cluster is not a valid cluster object.") 
    else if (!(require(snow))) 
      stop("Can't find required packages: snow")
    else if (!isClusterRunning(cluster))
      stop("the cluster is stopped.")
    else {

      # enter in cluster-aware mode
      cluster.aware = TRUE
      # set the test counter in all the cluster nodes
      clusterEvalQ(cluster, assign("test.counter", 0, envir = .GlobalEnv))
      # disable debugging, the slaves do not cat() here
      if (debug) {

        warning("disabling debugging output in cluster-aware mode.")
        debug = FALSE

      }#THEN

    }#ELSE

  }#THEN

  # NOTE: whitelist and blacklist relationship is the same as hosts.allow
  # and hosts.deny.

  if (!is.null(whitelist)) {

    if (class(whitelist) %in% c("matrix", "data.frame")) { 

       if (dim(whitelist)[2] != 2)
         stop("whitelist must have two columns.")

       if (is.data.frame(whitelist))
         whitelist = as.matrix(cbind(as.character(whitelist[,1]), 
           as.character(whitelist[,2])))

    }#THEN
    else if (is.character(whitelist)) {

      if (length(whitelist) != 2) 
        stop("whitelist must have two columns.") 

      whitelist = matrix(whitelist, ncol = 2, byrow = TRUE)

    }#THEN
    else {

      stop("whitelist must be a matrix or data.frame with two columns.")

    }#ELSE

    # add column names for easy reference.
    colnames(whitelist) = c("from", "to")
    # drop duplicate rows.
    unique(whitelist)

    # check all the names in the whitelist against the column names of x.
    if (any(!(unique(as.vector(whitelist)) %in% names(x))))
      stop("")

  }#THEN

  if (!is.null(blacklist)) {

    if (class(blacklist) %in% c("matrix", "data.frame")) { 

       if (dim(blacklist)[2] != 2) 
         stop("blacklist must have two columns.")

       if (is.data.frame(blacklist))
         blacklist = as.matrix(cbind(as.character(blacklist[,1]), 
           as.character(blacklist[,2])))

    }#THEN
    else if (is.character(blacklist)) {

      if (length(blacklist) != 2) 
        stop("blacklist must have two columns.") 

      blacklist = matrix(blacklist, ncol = 2, byrow = TRUE)

    }#THEN
    else {

      stop("blacklist must be a matrix or data.frame with two columns.")

    }#ELSE

    # add column names for easy reference.
    colnames(blacklist) = c("from", "to")
    # drop duplicate rows.
    unique(blacklist)

    # check all the names in the whitelist against the column names of x.
    if (any(!(unique(as.vector(blacklist)) %in% names(x))))
      stop("")

  }#THEN

  # if x -> y is whitelisted but y -> x is not, it is to be blacklisted.
  if (!is.null(whitelist)) {

    apply(whitelist, 1, 
      function(x) { 
        if (!is.whitelisted(whitelist, x[c(2,1)]))
          assign("blacklist", rbind(blacklist, x[c(2,1)]), 
            envir = sys.frame(-2))
      }#FUNCTION
    )

  }#THEN

  # if x -> y is whitelisted, it is to be removed from the blacklist.
  if (!is.null(blacklist) && !is.null(whitelist)) {

    blacklist = blacklist[!apply(blacklist, 1, 
      function(x){ is.whitelisted(whitelist, x) }),]

    blacklist = matrix(blacklist, ncol = 2, byrow = FALSE,
      dimnames = list(NULL, c("from", "to")))

  }#THEN

  # call the grow-shrink backend
  if (method == "gs") {

    if (cluster.aware) {
 
      res = grow.shrink.cluster(x = x, cluster = cluster, 
        whitelist = whitelist, blacklist = blacklist, test = test, 
        alpha = alpha, strict = strict, direction = direction, 
        debug = debug)

    }#THEN
    else if (optimized) {

      res = grow.shrink.optimized(x = x, whitelist = whitelist,
        blacklist = blacklist, test = test, alpha = alpha,
        strict = strict, direction = direction, debug = debug)

    }#THEN
    else {

      res = grow.shrink(x = x, whitelist = whitelist, blacklist = blacklist, 
        test = test, alpha = alpha, strict = strict, direction = direction, 
        debug = debug)

    }#ELSE

  }#THEN
  else if (method == "iamb") {

    if (cluster.aware) {
 
      res = incremental.association.cluster(x = x, cluster = cluster, 
        whitelist = whitelist, blacklist = blacklist, test = test, 
        alpha = alpha, strict = strict, direction = direction, 
        debug = debug)

    }#THEN
    else if (optimized) {

      res = incremental.association.optimized(x = x, whitelist = whitelist, 
        blacklist = blacklist, test = test, alpha = alpha, strict = strict,
        direction = direction, debug = debug)

    }#THEN
    else {

      res = incremental.association(x = x, whitelist = whitelist,
        blacklist = blacklist, test = test, alpha = alpha, strict = strict, 
        direction = direction, debug = debug)

    }#ELSE

  }#THEN
  else if (method == "fast-iamb") {

    if (cluster.aware) {

      res = fast.incremental.association.cluster(x = x, cluster = cluster,
        whitelist = whitelist, blacklist = blacklist, test = test,
        alpha = alpha, strict = strict, direction = direction,
        debug = debug)

    }#THEN
    else if (optimized) {

      res = fast.incremental.association.optimized(x = x, whitelist = whitelist,
        blacklist = blacklist, test = test, alpha = alpha, strict = strict, 
        direction = direction, debug = debug)

    }#THEN
    else {

      res = fast.incremental.association(x = x, whitelist = whitelist, 
        blacklist = blacklist, test = test, alpha = alpha, strict = strict, 
        direction = direction, debug = debug)

    }#ELSE

  }#THEN
  else if (method == "inter-iamb") {

    if (cluster.aware) {

      res = inter.incremental.association.cluster(x = x, cluster = cluster,
        whitelist = whitelist, blacklist = blacklist, test = test,
        alpha = alpha, strict = strict, direction = direction,
        debug = debug)

    }#THEN
    else if (optimized) {

      res = inter.incremental.association.optimized(x = x, whitelist = whitelist, 
        blacklist = blacklist, test = test, alpha = alpha, strict = strict, 
        direction = direction, debug = debug)

    }#THEN
    else {

      res = inter.incremental.association(x = x, whitelist = whitelist, 
        blacklist = blacklist, test = test, alpha = alpha, strict = strict, 
        direction = direction, debug = debug)

    }#ELSE

  }#THEN

  # add tests performed by the slaves to the test counter
  if (cluster.aware)
    res$ntests = res$ntests + 
      sum(unlist(clusterEvalQ(cluster, get("test.counter", envir = .GlobalEnv))))
  # save the learning method used.
  res$algo = method
  # save a discrete/continuous boolean value.
  res$discrete = test %in% available.discrete.tests
  # set the class of the return value.
  class(res) = "bn"

  invisible(res)

}#BNLEARN

