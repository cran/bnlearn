
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
  # check the number of levels of discrete variables, to guarantee that
  # the degrees of freedom of the tests are positive.
  if (is.data.discrete(x))
    for (col in names(x))
      if (nlevels(x[, col]) < 2)
        stop("all factors must have at leat two levels.")

}#CHECK.DATA

# check a single node label of a graph.
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

# check nodes (not from a bn object).
check.nodes = function(nodes, min.nodes = NULL) {

  # nodes must be a vector of character strings.
  if (!is(nodes, "character"))
    stop("nodes must be a vector of character strings, the labels of the nodes.")
  # at least one node is needed.
  if (length(nodes) < 1)
    stop("at leat one node label is needed.")
  # no duplicates allowed.
  if (any(duplicated(nodes)))
     stop("node labels must be unique.")
  # no empty strings
  if (any(nodes == ""))
    stop("an empty string is not a valid node label.")
  # minimum number of nodes requirement.
  if (!is.null(min.nodes))
    if (length(nodes) < min.nodes)
      stop(paste("at least", min.nodes, "nodes needed."))

}#CHECK.NODES

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

  # drop duplicate rows.
  whitelist = unique(whitelist)
  # add column names for easy reference.
  colnames(whitelist) = c("from", "to")

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

    # drop duplicate rows.
    blacklist = unique(blacklist)
    # add column names for easy reference.
    colnames(blacklist) = c("from", "to")

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
    if (!is.null(blacklist)) {

      blacklist = blacklist[!apply(blacklist, 1,
        function(x){ is.whitelisted(whitelist, x) }),]

      blacklist = matrix(blacklist, ncol = 2, byrow = FALSE,
        dimnames = list(NULL, c("from", "to")))

    }#THEN

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

  for (i in 1:ncol(data))
    if (!is(data[, i], "factor"))
      return(FALSE)

  return(TRUE)

}#IS.DATA.DISCRETE

# will the bayesian network be a continuous one?
is.data.continuous = function(data) {

  for (i in 1:ncol(data))
    if (!is(data[, i], "numeric"))
      return(FALSE)

  return(TRUE)

}#IS.DATA.CONTINUOUS

# there are missing data?
missing.data = function(data) {

  !all(complete.cases(data))

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

    # check if there is an imaginary sample size stored in the bn object;
    # otherwise use a the lowest possible value (3) if the network is
    # empty of the number of its parameters it is not.
    if (!is.null(network$learning$args$iss)) {

      iss = network$learning$args$iss

    }#THEN
    else if (nrow(network$arcs) == 0) {

      # even if the network is empty, the score should scale with the
      # number of nodes.
      iss = ncol(data)

    }#THEN
    else {

      # if the network is not empty, use a somewhat larger imaginary sample
      # size to give it and adequate (but still small) weight.
      if (is.data.discrete(data))
        iss = sum(nparams.backend(network, data, real = FALSE))
      else
        iss = ncol(data) * (ncol(data) + 1)/2

    }#THEN

  }#ELSE

return(iss)

}#CHECK.ISS

# check the phi defintion to be used in the bge score.
check.phi = function(phi, network, data) {

  if (!is.null(phi)) {

    if (!(phi %in% c("heckerman", "bottcher")))
      stop("unknown phi definition, should be either 'heckerman' or 'bottcher'.")

  }#THEN
  else {

    # check if there is an phi definition stored in the bn object;
    # otherwise use the one by heckerman.
    if (!is.null(network$learning$args$phi))
      phi = network$learning$args$phi
    else
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

# sanitize the extra arguments passed to the network scores.
check.score.args = function(score, network, data, extra.args) {

  if (score %in% c("dir", "bde")) {

    # check the imaginary sample size.
    extra.args$iss = check.iss(iss = extra.args$iss,
      network = network, data = data)

  }#THEN
  else if (score %in% c("aic", "bic")) {

    if (!is.null(extra.args$k)) {

      # validate the penalty weight.
      if (!is.positive(extra.args$k))
        stop("the penalty weight must be a positive numeric value.")

    }#THEN
    else {

      # set the penalty according to the chosen score.
      if (score == "aic") extra.args$k = 1
      else extra.args$k = log(nrow(data))/2

    }#ELSE

  }#THEN
  else if (score == "bge") {

    # check the imaginary sample size.
    extra.args$iss = check.iss(iss = extra.args$iss,
      network = network, data = data)

    # check phi estimator.
    extra.args$phi = check.phi(phi = extra.args$phi,
      network = network, data = data)

  }#THEN

  check.unused.args(extra.args, score.extra.args[[score]])

  return(extra.args)

}#CHECK.SCORE.ARGS

# warn about unused arguments.
check.unused.args = function(dots, used.args) {

  unused.args = !(names(dots) %in% used.args)
  if (any(unused.args))
    warning(paste("unused argument(s):", paste(names(dots)[unused.args],
              sep = "", collapse = " ")))

}#CHECK.UNUSED.ARGS

# check the the target nominal type I error rate
check.alpha = function(alpha, network = NULL) {

  # check the the target nominal type I error rate
  if (!is.null(alpha)) {

    # validate alpha.
    if (!is.positive(alpha) || (alpha > 1))
      stop("alpha must be a numerical value in [0,1].")

  }#THEN
  else {

    # check if there is an alpha value stored in the bn object;
    # otherwise use the usual 0.05 value.
    if (!is.null(network$learning$args$alpha))
      alpha = network$learning$args$alpha
    else
      alpha = 0.05

  }#ELSE

  return(alpha)

}#CHECK.ALPHA

check.amat = function(amat, nodes) {

  # a node is needed.
  if (missing(amat))
    stop("no adjacency matrix specified.")
  # the adjacency matrix must, well, be a matrix.
  if (!is(amat, "matrix"))
    stop("this is not a valid adjacency matrix.")
  # column names must be valid node labels.
  if (!is.null(colnames(amat)))
    if (!all(colnames(amat) %in% nodes))
      stop("node (column label) not present in the graph.")
  # column names must be valid node labels.
  if (!is.null(rownames(amat)))
    if (!all(rownames(amat) %in% nodes))
      stop("node (row label) not present in the graph.")
  # check the elements of the matrix.
  if (!all(amat %in% 0:1))
    stop("all the elements of an adjacency matrix must be equal to either 0 or 1.")
  # no arcs from a node to itself.
  if(any(diag(amat) != 0))
    stop("the elements on the diagonal must be zero.")

  return(amat)

}#CHECK.AMAT

# check logical flags.
check.logical = function(bool) {

  if (!is.logical(bool) || is.na(bool)) {

    stop(sprintf("%s must be a logical value (TRUE/FALSE).",
           deparse(substitute(bool))))

  }#THEN

}#CHECK.LOGICAL

# check logical flags.
check.bn = function(bn) {

  if (!is(bn, "bn")) {

    stop(sprintf("%s must be an object of class 'bn'.",
           deparse(substitute(bn))))

  }#THEN

}#CHECK.BN

# check parameters related to the random restart functions.
check.restart = function(restart, perturb) {

  if (!is.positive(restart)) {

    if ((restart != 0) || (length(restart) != 1))
      stop("the number of random restarts must be a non-negative numeric value.")

  }#THEN
  else if (!is.positive(perturb)) {

    stop("the number of changes at each radom restart must be a non-negative numeric value.")

  }#THEN

}#CHECK.RESTART

check.bn.vs.data = function(bn, data) {

  # the number of variables must be the same
  if (length(names(bn$nodes)) != ncol(data))
    stop("the network and the data have different numbers of variables.")
  # the variables must be the same.
  if (length(setdiff(names(bn$nodes) , names(data))) != 0)
    stop("the variables in the data and in the network do not match.")
  # data type versus network type
  if (bn$learning$test %in% c(available.discrete.tests, available.discrete.scores) &&
      is.data.continuous(data))
    stop("continuous data and discrete network.")
  if (bn$learning$test %in% c(available.continuous.tests, available.continuous.scores) &&
      is.data.discrete(data))
    stop("discrete data and continuous network.")

}#CHECK.BN.VS.DATA
