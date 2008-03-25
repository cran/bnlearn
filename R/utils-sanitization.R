
# is x a positive number?
is.positive = function(x) {

  is.numeric(x) &&
  (length(x) == 1) &&
  is.finite(x) &&
  (x > 0)

}#IS.POSITIVE

# check the data set.
check.data = function(x) {

  # check the data are there.
  if (missing(x))
    stop("the data are missing.")
  # x must be a data frame.
  if(!is.data.frame(x))
    stop("the data must be in a data frame.")
  # check the data for NULL/NaN/NA.
  if (missing.data(x))
    stop("the data set contains NULL/NaN/NA values.")
  # check the variables are either all continuous or all discrete.
  if (!is.data.discrete(x) && !is.data.continuous(x))
    stop("variables must be either all numeric or all factors.")

}#CHECK.DATA

# check a single node label.
check.node = function(node, graph) {

  # a node is needed.
  if (missing(node))
    stop("no node specified.")
  # a node label must be a character string.
  if (!is(node, "character"))
    stop("node must be a character string, the label of a node.")
  # only one node is needed.
  if (length(node) > 1)
    stop("only a single node may be specified.")
  # node must be a valid node label.
  if (!(node %in% names(graph$nodes)))
    stop("node not present in the graph.")

}#CHECK.NODE

# build a valid whitelist.
build.whitelist = function(whitelist, nodes) {

  if (is.null(whitelist)) {

    # no whitelist, nothing to do.
    return(NULL)

  }#THEN

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
  whitelist = unique(whitelist)

  # check all the names in the whitelist against the column names of x.
  if (any(!(unique(as.vector(whitelist)) %in% nodes)))
    stop("unknown node label present in the whitelist.")

  # if the whitelist itself contains cycles, no acyclic graph
  # can be learned.
  if (!is.acyclic(whitelist, nodes))
    stop("this whitelist does not allow an acyclic graph.")

whitelist

}#BUILD.WHITELIST

# build a valid blacklist.
build.blacklist = function(blacklist, whitelist, nodes) {

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
    blacklist = unique(blacklist)

    # check all the names in the blacklist against the column names of x.
    if (any(!(unique(as.vector(blacklist)) %in% nodes)))
      stop("unknown node label present in the blacklist.")

  }#THEN

  # update blacklist to agree with whitelist.
  # NOTE: whitelist and blacklist relationship is the same as hosts.allow
  # and hosts.deny.
  if (!is.null(whitelist)) {

    # if x -> y is whitelisted but y -> x is not, it is to be blacklisted.
    apply(whitelist, 1,
      function(x) {
        if (!is.whitelisted(whitelist, x[c(2,1)]))
          assign("blacklist", rbind(blacklist, x[c(2,1)]),
            envir = sys.frame(-2))
      })

    # if x -> y is whitelisted, it is to be removed from the blacklist.
    blacklist = blacklist[!apply(blacklist, 1,
      function(x){ is.whitelisted(whitelist, x) }),]

    blacklist = matrix(blacklist, ncol = 2, byrow = FALSE,
      dimnames = list(NULL, c("from", "to")))

  }#THEN

blacklist

}#BUILD.BLACKLIST

# check score labels.
check.score = function(score, data) {

  if (!is.null(score)) {

    # check the score/test label.
    if (!(score %in% available.scores))
      stop(paste("valid scores are:",
             paste(available.scores, collapse = " ")))
    # check if it's the right score for the data (discrete, continuous).
    if (!is.data.discrete(data) && (score %in% available.discrete.scores))
      stop(paste("score '", score, "' may be used with discrete data only.", sep = ""))
    if (!is.data.continuous(data) && (score %in% available.continuous.scores))
      stop(paste("score '", score, "' may be used with continuous data only.", sep = ""))

    return(score)

  }#THEN
  else {

    if (is.data.discrete(data))
      return("aic")
    else
      return("bge")

  }#ELSE

}#CHECK.SCORE

# check test labels.
check.test = function(test, data) {

  if (!is.null(test)) {

    # check the score/test label.
    if (!(test %in% available.tests))
      stop(paste("valid tests are:",
             paste(available.tests, collapse = " ")))
    # check if it's the right test for the data (discrete, continuous).
    if (!is.data.discrete(data) && (test %in% available.discrete.tests))
      stop(paste("test '", test, "' may be used with discrete data only.", sep = ""))
    if (!is.data.continuous(data) && (test %in% available.continuous.tests))
      stop(paste("test '", test, "' may be used with continuous data only.", sep = ""))

    return(test)

  }#THEN
  else {

    if (is.data.discrete(data))
      return("mi")
    else
      return("cor")

  }#ELSE

}#CHECK.TEST

# will the bayesian network be a discrete one?
is.data.discrete = function(data) {

  all(sapply(data, class) == "factor")

}#IS.DATA.DISCRETE

# will the bayesian network be a continuous one?
is.data.continuous = function(data) {

  all(sapply(data, class) == "numeric")

}#IS.DATA.CONTINUOUS

# there are missing data?
missing.data = function(data) {

  any(mapply(function(x) {is.na(x) || is.nan(x) || is.null(x)}, data))

}#MISSING.DATA

# check the imaginary sample size
check.iss = function(iss, network, data) {

  if (!is.null(iss)) {

    # validate the imaginary sample size.
    if (!is.positive(iss))
      stop("the imaginary sample size must be a positive numeric value.")
    # if iss = 1 the bge is NaN, if iss = 2 and phi = "heckerman" the
    # computation stops with the following error:
    # Error in solve.default(phi[A, A]) :
    #   Lapack routine dgesv: system is exactly singular
    if(is.data.continuous(data) && (iss < 3))
      stop("the imaginary sample size must be at least 3.")

  }#THEN
  else {

    # set an imaginary sample size to populate the a priori distribution.
    if (is.data.discrete(data)) {

      iss = 2 * prod(sapply(names(network$nodes),
        function(node) { nlevels(data[, node]) }))

    }#THEN
    else {

      iss = ncol(data) * (ncol(data) + 1)/2

    }#ELSE

  }#ELSE

as.integer(iss)

}#CHECK.ISS

# check the phi defintion to be used in the bge score.
check.phi = function(phi, network, data) {

  if (!is.null(phi)) {

    if (!(phi %in% c("heckerman", "bottcher")))
      stop("unknown phi definition, should be either 'heckerman' or 'bottcher'.")

  }#THEN
  else {

    phi = "heckerman"

  }#ELSE

phi

}#CHECK.PHI

check.arc = function(arc, network) {

  # check the arc is there.
  if (missing(arc))
    stop("the arc is missing.")
 # nodes must be a vector of character strings.
  if (!is(arc, "character"))
    stop("arc must be a vector of character strings, the labels of the nodes.")
  # exactly two nodes are needed.
  if (length(arc) != 2)
    stop("an arc if formed by exactly 2 nodes.")
  # no duplicates allowed.
  if(any(duplicated(arc)))
     stop("node labels in the arc must be unique.")
  # both elements of the arc must be valid node labels.
  if (!all(arc %in% names(network$nodes)))
    stop("node not present in the graph.")

}#CHECK.ARC
