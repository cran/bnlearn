
# is x a positive number?
is.positive = function(x) {

  is.numeric(x) &&
  (length(x) == 1) &&
  is.finite(x) &&
  (x > 0)

}#IS.POSITIVE

# is x a non-negative number?
is.non.negative = function(x) {

  is.numeric(x) &&
  (length(x) == 1) &&
  is.finite(x) &&
  (x >= 0)

}#IS.NON.NEGATIVE

# is x a positive integer?
is.positive.integer = function(x) {

  is.positive(x) && ((x %/% 1) == x)

}#IS.POSITIVE.INTEGER

# is x a vector of positive numbers?
is.positive.vector = function(x) {

  is.numeric(x) &&
  all(is.finite(x)) &&
  all(x > 0)

}#IS.POSITIVE.VECTOR

# is x a vector of non-negative numbers?
is.nonnegative.vector = function(x) {

  is.numeric(x) &&
  all(is.finite(x)) &&
  all(x >= 0)

}#IS.NONNEGATIVE.VECTOR

# is x a probability?
is.probability = function(x) {

  is.numeric(x) &&
  (length(x) == 1) &&
  is.finite(x) &&
  (x >= 0) &&
  (x <= 1)

}#IS.PROBABILITY

# is x a vector of probabilities?
is.probability.vector = function(x) {

  is.numeric(x) &&
  all(is.finite(x)) &&
  all(x >= 0) &&
  all(x <= 1) &&
  any(x > 0)

}#IS.PROBABILITY.VECTOR

# is x a single character string?
is.string = function(x) {

  is.character(x) &&
  (length(x) == 1) &&
  (x != "")

}#IS.STRING

is.ndmatrix = function(x) {

  is(x, c("table", "matrix", "array"))

}#IS.NDMATRIX

# is x a symmetric matrix?
is.symmetric = function(x) {

  .Call("is_symmetric",
        matrix = x,
        PACKAGE = "bnlearn")

}#IS.SYMMETRIC

is.cauchy.schwarz = function(x) {

  .Call("is_cauchy_schwarz",
        matrix = x,
        PACKAGE = "bnlearn")

}#Is.CAUCHY.SCHWARZ

# check the data set.
check.data = function(x, allow.mixed = FALSE) {

  # check the data are there.
  if (missing(x))
    stop("the data are missing.")
  # x must be a data frame.
  if(!is.data.frame(x))
    stop("the data must be in a data frame.")
  # check the data for NULL/NaN/NA.
  if (missing.data(x))
    stop("the data set contains NULL/NaN/NA values.")
  if (allow.mixed) {

    # check whether the variables are all factors or numeric.
    if (!is.data.mixed(x))
      stop("variables must be either numeric or factors.")

  }#THEN
  else {

    # check whether the variables are either all continuous or all discrete.
    if (!is.data.discrete(x) && !is.data.continuous(x))
      stop("variables must be either all real numbers or all factors.")

  }
  # check the number of levels of discrete variables, to guarantee that
  # the degrees of freedom of the tests are positive.
  if (is.data.discrete(x))
    for (col in names(x))
      if (nlevels(x[, col]) < 2)
        stop("all factors must have at least two levels.")

}#CHECK.DATA

# check nodes (not necessarily from a bn object).
check.nodes = function(nodes, graph = NULL, min.nodes = 1, max.nodes = Inf) {

  # a node is needed.
  if (missing(nodes))
    stop("no node specified.")
  # nodes must be a vector of character strings.
  if (!is(nodes, "character"))
    stop("nodes must be a vector of character strings, the labels of the nodes.")
  # no duplicates allowed.
  if (any(duplicated(nodes)))
     stop("node labels must be unique.")
  # no empty strings.
  if (any(is.na(nodes)) || any(nodes == ""))
    stop("an empty string is not a valid node label.")
  # maximum number of nodes requirement.
  if (length(nodes) > max.nodes)
    stop(paste("at most", max.nodes, "node(s) needed."))
  # minimum number of nodes requirement (usually 1).
  if (length(nodes) < min.nodes)
    stop(paste("at least", min.nodes, "node(s) needed."))
  # node must be a valid node label.
  if (!is.null(graph)) {

    if (is(graph, "bn")) {

      if (!all(nodes %in% names(graph$nodes)))
        stop(paste(c("node(s)", nodes[!(nodes %in% names(graph$nodes))],
               "not present in the graph."), collapse = " "))

    }#THEN
    else if (is(graph, "bn.fit")) {

      if (!all(nodes %in% names(graph)))
        stop(paste(c("node(s)", nodes[!(nodes %in% names(graph))],
               "not present in the graph."), collapse = " "))

    }#THEN
    else if (is.character(graph)) {

      if (!all(nodes %in% graph))
        stop(paste(c("node(s)", nodes[!(nodes %in% graph)],
               "not present in the graph."), collapse = " "))

    }#THEN

  }#THEN

}#CHECK.NODES

# check an arc set.
check.arcs = function(arcs, nodes) {

  # sanitize the set of arcs.
  if (is(arcs, "matrix") || is(arcs, "data.frame")) {

     if (dim(arcs)[2] != 2)
       stop("arc sets must have two columns.")

     if (is.data.frame(arcs))
       arcs = as.matrix(cbind(as.character(arcs[, 1]),
         as.character(arcs[, 2])))

     # be sure to set the column names.
     dimnames(arcs) = list(c(), c("from", "to"))

  }#THEN
  else if (is.character(arcs)) {

    # if there is an even number of labels fit them into a 2-column matrix.
    if ((length(arcs) %% 2) != 0)
      stop("arc sets must have two columns.")

    arcs = matrix(arcs, ncol = 2, byrow = TRUE,
              dimnames = list(c(), c("from", "to")))

  }#THEN
  else {

     stop("an arc set must be a matrix or data.frame with two columns.")

  }#ELSE

  # nodes must be valid node labels.
  if (!all(arcs %in% nodes))
    stop(paste(c("node(s)", unique(arcs[!(arcs %in% nodes)]),
           "not present in the graph."), collapse = " "))

  # remove duplicate arcs.
  arcs = unique.arcs(arcs, nodes, warn = TRUE)

  # check there are no loops among the arcs.
  loop = (arcs[, "from"] == arcs[, "to"])

  if (any(loop))
    stop(paste(c("invalid arcs that are actually loops:\n",
      paste("  ", arcs[loop, 1], "->", arcs[loop, 2], "\n"))))

  return(arcs)

}#CHECK.ARCS

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
      whitelist = as.matrix(cbind(as.character(whitelist[, 1]),
        as.character(whitelist[, 2])))

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
  whitelist = unique.arcs(whitelist, nodes, warn = TRUE)
  # add column names for easy reference.
  colnames(whitelist) = c("from", "to")

  # check all the names in the whitelist against the column names of x.
  if (any(!(unique(as.vector(whitelist)) %in% nodes)))
    stop("unknown node label present in the whitelist.")

  # if the whitelist itself contains cycles, no acyclic graph
  # can be learned.
  if (!is.acyclic(whitelist, nodes))
    stop("this whitelist does not allow an acyclic graph.")

  return(whitelist)

}#BUILD.WHITELIST

# build a valid blacklist.
build.blacklist = function(blacklist, whitelist, nodes) {

  if (!is.null(blacklist)) {

    if (class(blacklist) %in% c("matrix", "data.frame")) {

      if (dim(blacklist)[2] != 2)
        stop("blacklist must have two columns.")

      if (is.data.frame(blacklist))
        blacklist = as.matrix(cbind(as.character(blacklist[, 1]),
          as.character(blacklist[, 2])))

    }#THEN
    else if (is.character(blacklist)) {

      if (length(blacklist) != 2)
        stop("blacklist must have two columns.")

      blacklist = matrix(blacklist, ncol = 2, byrow = TRUE)

    }#THEN
    else {

      stop("blacklist must be a matrix or data.frame with two columns.")

    }#ELSE

    # check all the names in the blacklist against the column names of x.
    if (any(!(unique(as.vector(blacklist)) %in% nodes)))
      stop("unknown node label present in the blacklist.")

    # drop duplicate rows.
    blacklist = unique.arcs(blacklist, nodes)

  }#THEN

  # update blacklist to agree with whitelist.
  # NOTE: whitelist and blacklist relationship is the same as hosts.allow
  # and hosts.deny.
  if (!is.null(whitelist)) {

    # if x -> y is whitelisted but y -> x is not, it is to be blacklisted.
    apply(whitelist, 1,
      function(x) {
        if (!is.whitelisted(whitelist, x[c(2, 1)]))
          assign("blacklist", rbind(blacklist, x[c(2, 1)]),
            envir = sys.frame(-2))
      })

    # if x -> y is whitelisted, it is to be removed from the blacklist.
    if (!is.null(blacklist)) {

      blacklist = blacklist[!apply(blacklist, 1,
        function(x){ is.whitelisted(whitelist, x) }),]

      blacklist = matrix(blacklist, ncol = 2, byrow = FALSE,
        dimnames = list(NULL, c("from", "to")))

      # drop duplicate rows.
      blacklist = unique.arcs(blacklist, nodes)

    }#THEN

  }#THEN

  return(blacklist)

}#BUILD.BLACKLIST

# check the list of networks passed to custom.strength().
check.customlist = function(custom, nodes) {

  # check
  if (!is(custom, "list"))
    stop("networks must be a list of objects of class 'bn' or of arc sets.")
  if(!all(sapply(custom, function(x) { is(x, "bn") || is(x, "matrix") })))
    stop("x must be a list of objects of class 'bn' or of arc sets.")

  validate = function(custom, nodes) {

    if (is(custom, "bn")) {

      check.nodes(names(custom$nodes), graph = nodes, min.nodes = length(nodes),
        max.nodes = length(nodes))

    }
    else if (is(custom, "matrix")) {

      check.arcs(arcs = custom, nodes = nodes)

    }#THEN
    else {

      stop("x must be a list of objects of class 'bn' or of arc sets.")

    }

  return(TRUE)

  }#VALIDATE

  if (!all(sapply(custom, validate, nodes = nodes)))
    stop("x must be a list of objects of class 'bn' or of arc sets.")

}#CHECK.CUSTOMLIST

# check score labels.
check.score = function(score, data) {

  if (!is.null(score)) {

    # check it's a single character string.
    check.string(score)
    # check the score/test label.
    if (!(score %in% available.scores))
      stop(paste(c("valid scores are:\n",
             sprintf("    %-10s %s\n", names(score.labels), score.labels)), sep = ""))
    # check if it's the right score for the data (discrete, continuous).
    if (!is.data.discrete(data) && (score %in% available.discrete.scores))
      stop(paste("score '", score, "' may be used with discrete data only.", sep = ""))
    if (!is.data.continuous(data) && (score %in% available.continuous.scores))
      stop(paste("score '", score, "' may be used with continuous data only.", sep = ""))

    return(score)

  }#THEN
  else {

    if (is.data.discrete(data))
      return("bic")
    else
      return("bic-g")

  }#ELSE

}#CHECK.SCORE

# check test labels.
check.test = function(test, data) {

  if (!is.null(test)) {

    # check it's a single character string.
    check.string(test)
    # check the score/test label.
    if (!(test %in% available.tests))
      stop(paste(c("valid tests are:\n",
             sprintf("    %-10s %s\n", names(test.labels), test.labels)), sep = ""))
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

check.criterion = function(criterion, data) {

  # check it's a single character string.
  check.string(criterion)
  # check criterion's label.
  if (criterion %in% available.tests)
    criterion = check.test(criterion, data)
  else if (criterion %in% available.scores)
    criterion = check.score(criterion, data)
  else
    stop(paste(c("valid tests are:\n",
      sprintf("    %-10s %s\n", names(test.labels), test.labels),
      "  valid scores are:\n",
      sprintf("    %-10s %s\n", names(score.labels), score.labels)),
      sep = ""))

  return(criterion)

}#CHECK.CRITERION

# check loss functions' labels.
check.loss = function(loss, data) {

  if (!is.null(loss)) {

    # check it's a single character string.
    check.string(loss)
    # check the score/test label.
    if (!(loss %in% loss.functions))
      stop(paste(c("valid loss functions are:\n",
             sprintf("    %-10s %s\n", names(loss.labels), loss.labels)), sep = ""))
    if (!is.data.discrete(data) && (loss %in% discrete.loss.functions))
      stop(paste("loss function '", loss, "' may be used with discrete data only.", sep = ""))
    if (!is.data.continuous(data) && (loss %in% continuous.loss.functions))
      stop(paste("loss function '", loss, "' may be used with continuous data only.", sep = ""))

    return(loss)

  }#THEN
  else {

    if (is.data.discrete(data))
      return("logl")
    else
      return("logl-g")

  }#ELSE

}#CHECK.LOSS

# check the method used to fit the parameters of the network.
check.fitting.method = function(method, data) {

  if (!is.null(method)) {

    # check it's a single character string.
    check.string(method)
    # check the score/test label.
    if (!(method %in% available.fitting.methods))
      stop(paste(c("valid fitting methods are:\n",
             sprintf("    %-10s %s\n", names(fitting.labels), fitting.labels)), sep = ""))
    # bayesian parameter estimation is implemented only for discrete data.
    if (is.data.continuous(data) && (method == "bayes"))
      stop("Bayesian parameter estimation for Gaussian Bayesian networks is not implemented.")

    return(method)

  }#THEN
  else {

      return("mle")

  }#ELSE

}#CHECK.FITTING.METHOD

# check the method used to discretize the data.
check.discretization.method = function(method) {

  if (!is.null(method)) {

    # check it's a single character string.
    check.string(method)
    # check the score/test label.
    if (!(method %in% available.discretization.methods))
      stop(paste(c("valid discretization methods are:\n",
             sprintf("    %-10s %s\n", names(discretization.labels), discretization.labels)), sep = ""))

    return(method)

  }#THEN
  else {

      return("quantile")

  }#ELSE

}#CHECK.DISCRETIZATION.METHOD

# check the estimator for the mutual information.
check.mi.estimator = function(estimator, data) {

   if (!is.null(estimator)) {

    # check it's a single character string.
    check.string(estimator)
    # check the score/estimator label.
    if (!(estimator %in% available.mi))
      stop(paste(c("valid estimators are:\n",
             sprintf("    %-10s %s\n", names(mi.estimator.labels), mi.estimator.labels)), sep = ""))
    # check if it's the right estimator for the data (discrete, continuous).
    if (!is.data.discrete(data) && (estimator %in% available.discrete.mi))
      stop(paste("estimator '", estimator, "' may be used with discrete data only.", sep = ""))
    if (!is.data.continuous(data) && (estimator %in% available.continuous.mi))
      stop(paste("estimator '", estimator, "' may be used with continuous data only.", sep = ""))

    return(estimator)

  }#THEN
  else {

    if (is.data.discrete(data))
      return("mi")
    else
      return("mi-g")

  }#ELSE

}#CHECK.MI.ESTIMATOR

# is the fitted bayesian network a discrete one?
is.fitted.discrete = function(fitted) {

  for (i in 1:length(fitted))
    if (class(fitted[[i]]) != "bn.fit.dnode")
      return(FALSE)

  return(TRUE)

}#IS.FITTED.DISCRETE

# is the fitted bayesian network a continuous one?
is.fitted.continuous = function(fitted) {

  for (i in 1:length(fitted))
    if (class(fitted[[i]]) != "bn.fit.gnode")
      return(FALSE)

  return(TRUE)

}#IS.FITTED.CONTINUOUS

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
    if (!is.double(data[, i]))
      return(FALSE)

  return(TRUE)

}#IS.DATA.CONTINUOUS

# does the data include coth discrete and continuous variables?
is.data.mixed = function(data) {

  for (i in 1:ncol(data))
    if (!is.double(data[, i]) && !is(data[, i], "factor"))
      return(FALSE)

  return(TRUE)

}#IS.DATA.MIXED

# there are missing data?
missing.data = function(data) {

  !all(complete.cases(data))

}#MISSING.DATA

# check the imaginary sample size.
check.iss = function(iss, network, data) {

  if (!is.null(iss)) {

    # validate the imaginary sample size.
    if (!is.positive(iss) || (iss < 1))
      stop("the imaginary sample size must be a numeric value greater than 1.")
    # if iss = 1 the bge is NaN, if iss = 2 and phi = "heckerman" the
    # computation stops with the following error:
    # Error in solve.default(phi[A, A]) :
    #   Lapack routine dgesv: system is exactly singular
    if(is.data.continuous(data) && (iss < 3))
      stop("the imaginary sample size must be a numeric value greater than 3.")

  }#THEN
  else {

    # check whether there is an imaginary sample size stored in the bn object;
    # otherwise use a the de facto standard value of 10.
    if (!is.null(network$learning$args$iss))
      iss = network$learning$args$iss
    else
      iss = 10

  }#ELSE

  # coerce iss to integer.
  return(as.integer(iss))

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

  return(phi)

}#CHECK.PHI

# check the experimental data list.
check.experimental = function(exp, network, data) {

  if (!is.null(exp)) {

    if (!is.list(exp))
      stop("experimental data must be specified via a list of indexes.")
    if (!all(names(exp) %in% names(data)) || (length(names(exp)) == 0))
      stop("unkown variables specified in the experimental data list.")
    for (var in names(exp)) {

      if (!is.positive.vector(exp[[var]]))
        stop("indexes of experimental data must be positive integer numbers.")
      if (any(duplicated(exp[[var]])))
        stop("duplicated indexes for experimental data.")
      if (any(exp[[var]] > length(data[, var])))
        stop("out of bounds indexes for experimental data.")

      # just kill empty elements.
      if (length(exp[[var]]) == 0)
        exp[[var]] = NULL
      # also, convert evetything to integers to make things simpler at the
      # C level.
      exp[[var]] = as.integer(exp[[var]])

    }#FOR

  }#THEN
  else {

    # check whether there is a list stored in the bn object; if no experimental
    # data is specified, return an empty list (which is the same as using the
    # plain BDe score).
    if (!is.null(network$learning$args$exp))
      exp = network$learning$args$exp
    else
      exp = structure(vector(ncol(data), mode = "list"), names = names(data))

  }#ELSE

  return(exp)

}#CHECK.EXPERIMENTAL

# check the penalty used in AIC and BIC.
check.penalty = function(k, network, data, score) {

  if (!is.null(k)) {

    # validate the penalty weight.
    if (!is.positive(k))
      stop("the penalty weight must be a positive numeric value.")

  }#THEN
  else {

    # check whether there is a penaltization coefficient stored in the bn object,
    # use the default for the score function otherwise.
    if (!is.null(network$learning$args$k))
      k = network$learning$args$k
    else
      k = ifelse((score %in% c("aic", "aic-g")), 1, log(nrow(data))/2)

  }#ELSE

  return(k)

}#CHECK.PENALTY

# sanitize the extra arguments passed to the network scores.
check.score.args = function(score, network, data, extra.args) {

  if (score %in% c("bde", "bdes")) {

    # check the imaginary sample size.
    extra.args$iss = check.iss(iss = extra.args$iss,
      network = network, data = data)

  }#THEN
  else if (score == "mbde") {

    # check the imaginary sample size.
    extra.args$iss = check.iss(iss = extra.args$iss,
      network = network, data = data)

    # check the list of the experimental observations in the data set.
    extra.args$exp = check.experimental(exp = extra.args$exp,
      network = network, data = data)

  }#THEN
  else if (score %in% c("aic", "bic", "aic-g", "bic-g")) {

    extra.args$k = check.penalty(k = extra.args$k, network = network,
      data = data, score = score)

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

# sanitize the extra arguments passed to the random graph generation algorithms.
check.graph.generation.args = function(method, nodes, extra.args) {

  if (method == "ordered") {

    if (!is.null(extra.args$prob)) {

      # prob must be numeric.
      if (!is.probability(extra.args$prob))
        stop("the branching probability must be a numeric value in [0,1].")

    }#THEN
    else {

      # this default produces graphs with about the same number of
      # arcs as there are nodes.
      extra.args$prob = 2 / (length(nodes) - 1)

    }#ELSE

  }#THEN
  else if (method %in% c("ic-dag", "melancon")) {

    if (!is.null(extra.args$every)) {

      if (!is.positive(extra.args$every))
        stop("'every' must be a positive integer number.")

    }#THEN
    else {

      extra.args$every = 1

    }#ELSE

    if (!is.null(extra.args$burn.in)) {

      if (!is.positive(extra.args$burn.in))
        stop("the burn in length must be a positive integer number.")

    }#THEN
    else {

      extra.args$burn.in = 6 * length(nodes)^2

    }#ELSE

    if (!is.null(extra.args$max.in.degree)) {

      if (!is.positive.integer(extra.args$max.in.degree))
        stop("the maximum in-degree must be a positive integer number.")

      if (extra.args$max.in.degree >= length(nodes)) {

        warning("a node cannot have an in-degree greater or equal to the number of nodes in the graph.")
        warning("the condition on the in-degree will be ignored.")

      }#THEN

    }#THEN
    else {

      extra.args$max.in.degree = Inf

    }#ELSE

    if (!is.null(extra.args$max.out.degree)) {

      if (!is.positive.integer(extra.args$max.out.degree))
        stop("the maximum out-degree must be a positive integer number.")

      if (extra.args$max.out.degree >= length(nodes)) {

        warning("a node cannot have an out-degree greater or equal to the number of nodes in the graph.")
        warning("the condition on the out-degree will be ignored.")

      }#THEN

    }#THEN
    else {

      extra.args$max.out.degree = Inf

    }#ELSE

    if (!is.null(extra.args$max.degree)) {

      if (!is.positive.integer(extra.args$max.degree))
        stop("the maximum out-degree must be a positive integer number.")

      if (is.finite(extra.args$max.in.degree) &&
          extra.args$max.in.degree > extra.args$max.degree)
        stop("the maximun in-degree must be lesser or equal to the maximum degree.")

      if (is.finite(extra.args$max.out.degree) &&
          extra.args$max.out.degree > extra.args$max.degree)
        stop("the maximun out-degree must be lesser or equal to the maximum degree.")

      if (extra.args$max.degree >= length(nodes)) {

        warning("a node cannot have a degree greater or equal to the number of nodes in the graph.")
        warning("the condition on the degree will be ignored.")

      }#THEN

    }#THEN
    else {

      extra.args$max.degree = Inf

    }#ELSE

  }#THEN

  check.unused.args(extra.args, graph.generation.extra.args[[method]])

  return(extra.args)

}#CHECK.GRAPH.GENERATION.ARGS

# check bootstrap arguments (when they are passed as variable length args).
check.bootstrap.args = function(extra.args, network, data) {

  # check the number of bootstrap replicates.
  extra.args$R = check.replicates(extra.args$R)
  # check the size of each bootstrap sample.
  extra.args$m = check.bootsize(extra.args$m, data)
  # check the learning algorithm.
  algorithm = check.learning.algorithm(extra.args[["algorithm"]], bn = network)
  # check the extra arguments for the learning algorithm.
  algorithm.args = check.learning.algorithm.args(extra.args[["algorithm.args"]],
                     algorithm = algorithm, bn = network)

  extra.args[["algorithm"]] = algorithm
  extra.args[["algorithm.args"]] = algorithm.args

  # remap additional arguments used in hybrid algorithms.
  if (algorithm %in% hybrid.algorithms) {

    # there's no need to sanitize these parameters, it's done either in
    # bnlearn() or in greedy.search() already.
    if (is.null(extra.args[["algorithm.args"]]$restrict))
      extra.args[["algorithm.args"]]$restrict = network$learning$restrict
    if (is.null(extra.args[["algorithm.args"]]$maximize))
      extra.args[["algorithm.args"]]$maximize = network$learning$maximize
    if (is.null(extra.args[["algorithm.args"]]$test))
      extra.args[["algorithm.args"]]$test = network$learning$rstest
    if (is.null(extra.args[["algorithm.args"]]$score))
      extra.args[["algorithm.args"]]$score = network$learning$maxscore

  }#THEN

  # warn about unused arguments.
  check.unused.args(extra.args, c("R", "m", "algorithm", "algorithm.args"))

  return(extra.args)

}#CHECK.BOOTSTRAP.ARGS

# sanitize the extra arguments passed to the conditional probability algorithms.
check.cpq.args = function(fitted, extra.args, method) {

  if (method %in% c("ls")) {

    if (!is.null(extra.args$n)) {

      if (!is.positive.integer(extra.args$n))
        stop("the number of observations to be sampled must be a positive integer number.")

    }#THEN
    else {

      # this is a rule of thumb, the error of the estimate has no closed-form
      # expression (Friedman & Koller).
      extra.args$n = 5000 * nparams.fitted(fitted)

    }#ELSE

    if (!is.null(extra.args$batch)) {

      if (!is.positive.integer(extra.args$batch))
        stop("the number of observations to be sampled must be a positive integer number.")

      if (extra.args$batch > extra.args$n) {

        warning("cannot generate a batch bigger than the whole generated data set.")
        warning("batch size will be ignored.")

      }#THEN

    }#THEN
    else {

      # perform small simulations in a single batch, and split larger ones.
      extra.args$batch = min(extra.args$n, 10^4)

    }#ELSE

  }#THEN

  # warn about unused arguments.
  check.unused.args(extra.args, cpq.extra.args[[method]])

  return(extra.args)

}#CHECK.CPQ.ARGS

# sanitize the extra arguments passed to loss functions.
check.loss.args = function(loss, bn, nodes, data, extra.args) {

  valid.args = loss.extra.args[[loss]]

  if (loss == "pred") {

    if (!is.null(extra.args$target)) {

      if (!is.string(extra.args$target) || !(extra.args$target %in% nodes))
        stop("target node must be a single, valid node label for the network.")

    }#THEN
    else {

      # the target node is obvious for classifiers.
      if (is(bn, c("bn.naive", "bn.tan"))) {

        if (is(bn, "bn"))
          extra.args$target = bn$learning$args$training
        else
          extra.args$target = attr(bn, "training")

      }#THEN
      else {

        stop("missing target node for which to compute the prediction error.")

      }#ELSE

    }#ELSE

    # check the prior distribution.
    if (is(bn, c("bn.naive", "bn.tan"))) {

      extra.args$prior = check.prior(extra.args$prior, data[, extra.args$target])
      valid.args = c(valid.args, "prior")

    }#THEN

  }#THEN

  # warn about unused arguments.
  check.unused.args(extra.args, valid.args)

  return(extra.args)

}#CHECK.LOSS.ARGS

# sanitize the extra arguments passed to fitting functions.
check.fitting.args = function(method, network, data, extra.args) {

  if (method == "bayes") {

    # check the imaginary sample size.
    extra.args$iss = check.iss(iss = extra.args$iss,
      network = network, data = data)

  }#THEN

  # warn about unused arguments.
  check.unused.args(extra.args, fitting.extra.args[[method]])

  return(extra.args)

}#CHECK.FITTING.ARGS

# sanitize the extra arguments passed to discretization methods.
check.discretization.args = function(method, data, breaks, extra.args) {

  if (method == "hartemink") {

    if (is.data.discrete(data)) {

      extra.args$ibreaks = nlevels(data[, 1])
      warning("data are already discrete, 'ibreaks' and 'idisc' are ignored.")

    }#THEN
    else {

      if (!is.null(extra.args$idisc)) {

        idisc = extra.args$idisc

        # check it's a single character string.
        check.string(idisc)
        # check the score/test label.
        other.methods = available.discretization.methods[available.discretization.methods != "hartemink"]
        if (!(idisc %in% other.methods))
          stop(paste(c("valid initial discretization methods are:\n",
                sprintf("    %-10s %s\n", other.methods, discretization.labels[other.methods])), sep = ""))


      }#THEN
      else {

        # default to quantile discretization as for Hartemink's recommendation.
        extra.args$idisc = "quantile"

      }#ELSE

      if (!is.null(extra.args$ibreaks)) {

        if (!is.positive.integer(extra.args$ibreaks))
          stop("the number of initial breaks must be a positive integer number.")
        if (extra.args$ibreaks < breaks)
          stop("insufficient number of levels, at least ", breaks, " required.")

      }#THEN
      else {

        ndata = nrow(data)

        if (ndata > 500)
          extra.args$ibreaks = 100
        else if (ndata > 100)
          extra.args$ibreaks = 50
        else if (ndata > 50)
          extra.args$ibreaks = 20
        else if (ndata > 10)
          extra.args$ibreaks = 10
        else
          extra.args$ibreaks = ndata

      }#ELSE

    }#ELSE

  }#THEN

  # warn about unused arguments.
  check.unused.args(extra.args, discretization.extra.args[[method]])

  return(extra.args)

}#CHECK.DISCRETIZATION.ARGS

# sanitize the extra arguments passed to Bayesian classifiers.
check.classifier.args = function(method, data, training, explanatory,
    extra.args) {

  if (method == "tan") {

    # check the label of the mutual information estimator.
    extra.args$estimator = check.mi.estimator(extra.args$estimator, data)

    # check the node to use the root of the tree (if not specified pick the first
    # explanatory variable assuming natural ordering).
    if (!is.null(extra.args$root))
      check.nodes(extra.args$root, graph = explanatory, max.nodes = 1)
    else
      extra.args$root = explanatory[1]

  }#THEN

  return(extra.args)

}#CHECK.CLASSIFICATION.ARGS

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
    if (!is.probability(alpha))
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

# check the number of permutation/boostrap samples.
check.B = function(B, criterion) {

  if (criterion %in% resampling.tests) {

    if (!is.null(B)) {

      if (!is.positive.integer(B))
        stop("the number of permutations/bootstrap replications must be a positive integer number.")

      B = as.integer(B)

    }#THEN
    else {

      B = 5000L

    }#ELSE

  }#THEN
  else {

    if (!is.null(B))
      warning("this test does not require any permutations/bootstrap resampling, ignoring B.\n")

    B = NULL

  }#ELSE

  return(B)

}#CHECK.B

check.amat = function(amat, nodes) {

  # a node is needed.
  if (missing(amat))
    stop("no adjacency matrix specified.")
  # the adjacency matrix must, well, be a matrix.
  if (!is(amat, "matrix") || (ncol(amat) != nrow(amat)) || (length(dim(amat)) != 2))
    stop("an adjacency matrix must be a 2-dimensional square matrix.")
  # check the dimensions against the number of nodes in the graph.
  if (any(dim(amat) != length(nodes)))
    stop("the dimensions of the adjacency matrix do not agree with the number of nodes in the graph.")
  # column names must be valid node labels.
  if (!is.null(colnames(amat)))
    if (!all(colnames(amat) %in% nodes))
      stop("node (column label) not present in the graph.")
  # column names must be valid node labels.
  if (!is.null(rownames(amat)))
    if (!all(rownames(amat) %in% nodes))
      stop("node (row label) not present in the graph.")
  # column names must match with row names.
  if (!is.null(colnames(amat)) && !is.null(rownames(amat))) {

    if (!identical(colnames(amat), rownames(amat)))
      stop("row/column names mismatch in the adjacency matrix.")

    if (!identical(colnames(amat), nodes) || !identical(rownames(amat), nodes)) {

      warning("rearranging the rows/columns of the adjacency matrix.")

      amat = amat[nodes, nodes, drop = FALSE]

    }#THEN

  }#THEN
  # make really sure the adjacency matrix is made up of integers.
  if (storage.mode(amat) != "integer")
    storage.mode(amat) = "integer"
  # check the elements of the matrix.
  if (!all((amat == 0L) | (amat == 1L)))
    stop("all the elements of an adjacency matrix must be equal to either 0 or 1.")
  # no arcs from a node to itself.
  if(any(diag(amat) != 0))
    stop("the elements on the diagonal must be zero.")

  return(amat)

}#CHECK.AMAT

check.covariance = function(m) {

  # the adjacency matrix must, well, be a matrix.
  if (!is(m, "matrix") || (ncol(m) != nrow(m)) || (length(dim(m)) != 2))
    stop("a covariance matrix must be a 2-dimensional square matrix.")
  # check the elements of the matrix.
  if (!is.numeric(m))
    stop("the elements of a covariance matrix must be real numbres.")
  # check whether the matrix is symmetric.
  if (!is.symmetric(m))
    stop("a covariance matrix must be symmetric.")
  # check whether the matrix obeys the Cauchy-Schwarz theorem.
  if (!is.cauchy.schwarz(m))
    stop("a covariance matrix must obey the Cauchy-Schwarz theorem.")

}#CHECK.COVARIANCE

# check logical flags.
check.logical = function(bool) {

  if (!is.logical(bool) || is.na(bool) || (length(bool) != 1)) {

    stop(sprintf("%s must be a logical value (TRUE/FALSE).",
           deparse(substitute(bool))))

  }#THEN

}#CHECK.LOGICAL

# check character strings.
check.string = function(string) {

  if (!is.string(string)) {

    stop(sprintf("%s must be a character string.",
           deparse(substitute(string))))

  }#THEN

}#CHECK.STRING

# check an object of class bn.
check.bn = function(bn) {

  if (missing(bn))
    stop("an object of class 'bn' is required.")
  if (!is(bn, "bn")) {

    stop(sprintf("%s must be an object of class 'bn'.",
           deparse(substitute(bn))))

  }#THEN

}#CHECK.BN

# check two bn's against each other.
match.bn = function(bn1, bn2) {

  # the two networks must have the same node set.
  nodes1 = names(bn1$nodes)
  nodes2 = names(bn2$nodes)

  equal = setequal(nodes1, nodes2) && (length(nodes1) == length(nodes2))

  if (!equal)
    stop("the two networks have different node sets.")

}#MATCH.BN

# check an object of class bn or bn.fit.
check.bn.or.fit = function(bn) {

  if (missing(bn))
    stop("an object of class 'bn' or 'bn.fit' is required.")
  if (!is(bn, "bn") && !is(bn, "bn.fit")) {

    stop(sprintf("%s must be an object of class 'bn' or 'bn.fit'.",
           deparse(substitute(bn))))

  }#THEN

}#CHECK.BN.OR.FIT

# check an object of class bn.
check.fit = function(bn) {

  if (missing(bn))
    stop("an object of class 'bn.fit' is required.")
  if (!is(bn, "bn.fit")) {

    stop(sprintf("%s must be an object of class 'bn.fit'.",
           deparse(substitute(bn))))

  }#THEN

}#CHECK.FIT

# check the structure of a naive Bayes classifier.
check.bn.naive = function(bn) {

  # check whether it's a valid bn/bn.fit object.
  check.bn.or.fit(bn)
  # there must be a single root node, check.
  root = root.leaf.nodes(bn, leaf = FALSE)

  if (length(root) != 1)
    stop("a naive Bayes classifier can have only one root node, the training variable.")

  # cache the node labels.
  if (is(bn, "bn"))
    nodes = names(bn$nodes)
  else
    nodes = names(bn)
  # get the explanatory variables.
  explanatory = nodes[nodes != root]
  leafs = root.leaf.nodes(bn, leaf = TRUE)
  # all the explanatory variables must be leaf nodes, check.
  if (!identical(sort(explanatory), sort(leafs)))
    stop("all the explanatory variables must be leaf nodes.")
  # all the explanatory variables must have a single parent, the root node, check.
  nparents = sapply(explanatory, function(node) { length(parents(bn, node))  })

  if (any(nparents != 1))
    stop("all the explanatory variables must be children of the training variable.")

}#CHECK.BN.NAIVE

# check the structure of a naive Bayes classifier.
check.bn.tan = function(bn) {

  # check whether it's a valid bn/bn.fit object.
  check.bn.or.fit(bn)
  # there must be a single root node, check.
  root = root.leaf.nodes(bn, leaf = FALSE)

  if (length(root) != 1)
    stop("a naive Bayes classifier can have only one root node, the training variable.")

  # that root node must be the training variable, check.
  if (is(bn, "bn")) {

    # double check just in case.
    check.nodes(bn$learning$args$training)

    nodes = names(bn$nodes)
    training = bn$learning$args$training

  }#THEN
  else {

    # double check just in case.
    check.nodes(attr(bn, "training"))

    nodes = names(bn)
    training = attr(bn, "training")

  }#ELSE

  if (!identical(training, root))
    stop("the training node is not the only root node in the graph.")

  # get the explanatory variables.
  explanatory = nodes[nodes != root]
  # all the explanatory variables save one must have exactly two parents, check.
  nparents = sapply(explanatory, function(node) { length(parents(bn, node))  })

  if (!( (length(which(nparents == 2)) == length(explanatory) - 1) && (length(which(nparents == 1)) == 1) ))
    stop("the explanatory variables must form a tree.")

}#CHECK.BN.TAN

# check an object of class bn.strength.
check.bn.strength = function(strength, nodes) {

  if (missing(strength))
    stop("an object of class 'bn.strength' is required.")
  if (!is(strength, "bn.strength")) {

    stop(sprintf("%s must be an object of class 'bn.strength'.",
           deparse(substitute(strength))))

  }#THEN
  if (!(ncol(strength) %in% 3:4))
    stop("objects of class 'bn.strength' must have 3 or 4 columns.")
  if (!identical(names(strength), c("from", "to", "strength")) &&
      !identical(names(strength), c("from", "to", "strength", "direction")))
    stop("objects of class 'bn.strength' must be data frames with column names ",
         "'from', 'to', 'strength' and (optionally) 'direction'.")
  if (!all(c("mode", "threshold") %in% names(attributes(strength))))
    stop("objects of class 'bn.strength' must have a 'mode' and a 'strength' attribute.")
  if (!missing(nodes))
    check.arcs(strength[, c("from", "to"), drop = FALSE], nodes)

}#CHECK.BN.STRENGTH

# sanitize the threshold value.
check.threshold = function(threshold, strength) {

  if (missing(threshold))
    threshold = attr(strength, "threshold")
  else {

    s = strength[, "strength"]

    if (!is.numeric(threshold) || (length(threshold) != 1) || is.nan(threshold))
      stop("the threshold must be a numeric value.")
    if ((threshold < min(s)) || (threshold > max(s)))
      warning("the threshold is outside the range of the strength values.")

  }#ELSE

  return(threshold)

}#CHECK.THRESHOLD

# check parameters related to the random restart functions.
check.restart = function(restart) {

  # set the default value if not specified.
  if (is.null(restart) || (restart == 0))
      return(0)

  if (!is.positive.integer(restart))
    stop("the number of random restarts must be a non-negative numeric value.")
  else
    return(restart)

}#CHECK.RESTART

check.perturb = function(perturb) {

  # set the default value if not specified.
  if (is.null(perturb))
      return(1)

  if (!is.positive.integer(perturb))
    stop("the number of changes at each radom restart must be a non-negative numeric value.")
  else
    return(perturb)

}#CHECK.PERTURB

# check the maximum number of iterations.
check.max.iter = function(max.iter) {

  # set the default value if not specified.
  if (is.null(max.iter))
    return(Inf)

  if ((max.iter != Inf) && !is.positive.integer(max.iter))
    stop("the maximum number of iterations must be a positive integer number.")
  else
    return(max.iter)

}#CHECK.MAX.ITER

# check arguments related to the tabu list.
check.tabu = function(tabu) {

  # set the default value if not specified.
  if (is.null(tabu))
    return(10)

  if (!is.positive.integer(tabu))
    stop("the length of the tabu list must be a positive integer number.")
  else
    return(tabu)

}#CHECK.TABU

check.max.tabu = function(max, tabu) {

  if (is.null(max))
    return(tabu)

  # check the number of iterations the algorithm can perform without
  # improving the best network score.
  if (!is.positive.integer(max))
    stop("the maximum number of iterations without any score improvement must be a positive integer number.")
  # the tabu list should be longer than that, otherwise the search can do a
  # U-turn and return to the local maximum it left before (thus creating an
  # endless loop).
  if (max > tabu)
    stop("the maximum number of iterations without any score improvement must not be grater than the length of the tabu list.")

  return(max)

}#CHECK.MAX.TABU

# check bn metadata against the data it's used with.
check.bn.vs.data = function(bn, data) {

  # the number of variables must be the same
  if (length(names(bn$nodes)) != ncol(data))
    stop("the network and the data have different numbers of variables.")
  # the variables must be the same.
  if (length(setdiff(names(bn$nodes) , names(data))) != 0)
    stop("the variables in the data and in the network do not match.")
  # data type versus network type.
  if (bn$learning$test %in% c(available.discrete.tests, available.discrete.scores) &&
      is.data.continuous(data))
    stop("continuous data and discrete network.")
  if (bn$learning$test %in% c(available.continuous.tests, available.continuous.scores) &&
      is.data.discrete(data))
    stop("discrete data and continuous network.")

}#CHECK.BN.VS.DATA

# check bn.fit metadata against the data it's used with.
check.fit.vs.data = function(fitted, data) {

  # the number of variables must be the same
  if (length(names(fitted)) != ncol(data))
    stop("the network and the data have different numbers of variables.")
  # the variables must be the same.
  if (length(setdiff(names(fitted) , names(data))) != 0)
    stop("the variables in the data and in the network do not match.")
  # data type versus network type.
  if (is.fitted.discrete(fitted)) {

    if (is.data.continuous(data))
      stop("continuous data and discrete network.")

    # double-check the levels of the variables against those of the nodes.
    for (node in names(fitted)) {

      data.levels = levels(data[, node])
      node.levels = dimnames(fitted[[node]]$prob)[[1]]

      if(!identical(data.levels, node.levels))
        stop("the levels of node '", node, "' do not match the levels of the ",
             "corresponding variable in the data.")

    }#FOR

  }#THEN
  if (is.fitted.continuous(fitted) && is.data.discrete(data))
    stop("discrete data and continuous network.")

}#CHECK.FIT.VS.DATA

# check bn.fit.{d,g}node metadata against the data it's used with.
check.fit.node.vs.data = function(fitted, data) {

  relevant = c(fitted$node, fitted$parents)

  # check whether all relevant nodes are in the data.
  if (!all(relevant %in% names(data)))
    stop("not all required nodes are present in the data.")
  # data type versus network type.
  if ((class(fitted) == "bn.fit.dnode") && is.data.continuous(data))
      stop("continuous data and discrete network.")
  if ((class(fitted) == "bn.fit.gnode") && is.data.discrete(data))
    stop("discrete data and continuous network.")
  # double-check the levels of the variables against those of the nodes.
  if (class(fitted) == "bn.fit.dnode") {

  for (node in relevant) {

      data.levels = levels(data[, node])
      node.levels = dimnames(fitted$prob)[[node]]

      if(!identical(data.levels, node.levels))
        stop("the levels of node '", node, "' do not match the levels of the ",
             "corresponding variable in the data.")

    }#FOR

  }#THEN

}#CHECK.FIT.NODE.VS.DATA

# check a colour identifier (not necessarily a string/integer).
check.colour = function(col) {

  if (identical(tryCatch(col2rgb(col), error = function(x) { FALSE }), FALSE))
    stop(sprintf("%s is not a valid colour identifier.",
           deparse(substitute(col))))

}#CHECK.COLOUR

# check the line type identifier.
check.lty = function(lty) {

  lty.strings = c("blank", "solid", "dashed", "dotted", "dotdash", "longdash",
                  "twodash")

  if (!(lty %in% 0:6) && !(lty %in% lty.strings))
    stop(sprintf("%s is not a valid line type identifier.",
           deparse(substitute(lty))))

}#CHECK.LTY

# check the label of a learning algorithm.
check.learning.algorithm = function(algorithm, class = "all", bn) {

  ok = character(0)

  if (missing(algorithm) || is.null(algorithm)) {

    # use the one specified by the bn object as the default.
    if (missing(bn))
      stop("the learning algorithm must be a character string.")
    else if (is(bn, "bn"))
      algorithm = bn$learning$algo

  }#THEN
  else if (!is.string(algorithm))
    stop("the learning algorithm must be a character string.")

  # select the right class of algorithms.
  if ("constraint" %in% class)
    ok = c(ok, constraint.based.algorithms)
  if ("markov.blanket" %in% class)
    ok = c(ok, markov.blanket.algorithms)
  if ("neighbours" %in% class)
    ok = c(ok, local.search.algorithms)
  if ("score" %in% class)
    ok = c(ok, score.based.algorithms)
  if ("mim" %in% class)
    ok = c(ok, mim.based.algorithms)
  if ("classifier" %in% class)
    ok = c(ok, classifiers)
  if ("all" %in% class)
    ok = available.learning.algorithms

  if (!(algorithm %in% ok))
       stop(paste(c("valid learning algorithms are:\n",
            sprintf("    %-10s %s\n", ok, method.labels[ok])), sep = ""))

  return(algorithm)

}#CHECK.LEARNING.ALGORITHM

# check the aruments of a learning algorithm (for use in bootstrap).
check.learning.algorithm.args = function(args, algorithm, bn) {

  # convert args into a list, if it's not one already.
  if (!is.list(args))
      args = as.list(args)

  # if a reference bn is specified, guess as many parameters as possbile.
  if (!(missing(algorithm) || missing(bn))) {

    # use the same score/conditional independence test.
    if (algorithm %in% constraint.based.algorithms) {

      # it's essential to check it's actually an independence test,
      # it could be a score function or NA.
      if (!("test" %in% names(args)))
        if (bn$learning$test %in% available.tests)
          args$test = bn$learning$test

      # set the appropriate value for the optimization flag.
      if (!("optimized" %in% names(args)))
        args$optimized = bn$learning$optimized

      # pass along all the parameters in bn$learning$args.
      if (length(bn$learning$args) > 0) {

        if (!("alpha" %in% names(args)))
          args$alpha = bn$learning$args$alpha

        if ("test" %in% names(args) && !("B" %in% names(args)))
          if (args$test %in% resampling.tests)
            args$B = bn$learning$args$B

      }#THEN

    }#THEN
    else if (algorithm %in% score.based.algorithms) {

      if (!("score" %in% names(args)))
        if (bn$learning$test %in% available.scores)
          args$score = bn$learning$test

      # set the appropriate value for the optimization flag.
      if (!("optimized" %in% names(args)))
        args$optimized = bn$learning$optimized

      # pass along the relevant parameters in bn$learning$args if the score
      # function is the same (hint: different scores have paramenters with
      # the same name but different meanings).
      if (("score" %in% names(args)) && (args$score == bn$learning$test))
        for (arg in names(bn$learning$args))
          if (!(arg %in% names(args)) && (arg %in% (score.extra.args[[args$score]])))
            args[[arg]] = bn$learning$args[[arg]]

    }#THEN

    # pass along whitelist and blacklist.
    if (!is.null(bn$learning$whitelist))
      args$whitelist = bn$learning$whitelist
    if (!is.null(bn$learning$blacklist))
      args$blacklist = bn$learning$blacklist

  }#THEN

  # remove any spurious x arguments, the data are provided by the bootstrap.
  if ("x" %in% names(args)) {

    args$x = NULL

    warning("removing 'x' from 'algorithm.args', the data set is provided by the bootstrap sampling.")

  }#THEN

  return(args)

}#CHECK.LEARNING.ALGORITHM.ARGS

# check the number of bootstrap replicates.
check.replicates = function(R, default = 200) {

  if (missing(R) || is.null(R))
    R = default
  else if (!is.positive.integer(R))
    stop("the number of bootstrap replicates must be a positive integer.")

  return(R)

}#CHECK.RESAMPLING

# check the size of bootstrap replicates.
check.bootsize = function(m, data, default = nrow(data)) {

  if (missing(m) || is.null(m))
    m = default
  else if (!is.positive.integer(m))
    stop("bootstrap sample size must be a positive integer.")

  return(m)

}#CHECK.BOOTSIZE

# check the label of the multivariate Bernulli variance test.
check.mvber.vartest = function(method) {

  if (missing(method))
    stop(paste(c("valid statistical tests are:\n",
           sprintf("    %-10s %s\n", names(mvber.labels), mvber.labels)), sep = ""))

  if (method %in% available.mvber.vartests)
    method = which(available.mvber.vartests %in% method)
  else
    stop(paste(c("valid statistical tests are:\n",
           sprintf("    %-10s %s\n", names(mvber.labels), mvber.labels)), sep = ""))

  return(method)

}#CHECK.MVBER.VARTEST

# check a prior distribution against the observed variable.
check.prior = function(prior, training) {

  if (is.factor(training)) {

    if (missing(prior) || is.null(prior)) {

      # use an empirical prior if none is provided by the user.
      prior = summary(training)

    }#THEN
    else {

      if (length(prior) != nlevels(training))
        stop("the prior distribution and the training variable have a different number of levels.")
      if (!is.probability.vector(prior))
        stop("the prior distribution must be expressed as a probability vector.")

    }#ELSE

  }#THEN

  # make sure the prior probabilities sum to one.
  prior = prior / sum(prior)

  return(prior)

}#CHECK.PRIOR

# check a vector of weights.
check.weights = function(weights, len) {

  if (missing(weights) || is.null(weights)) {

    weights = rep(1, len)

  }#THEN
  else {

    if (!is.nonnegative.vector(weights))
      stop("missing or negative weights are not allowed.")

    if (length(weights) != len)
      stop("wrong number of weights, ", length(weights),
        " weights while ", len, " are needed.")

    weights = prop.table(weights)

  }#ELSE

  return(weights)

}#CHECK.WEIGHTS

# check a user-specified bn.fit.gnode.
check.fit.gnode.spec = function(x) {

  components =  c("coef", "fitted", "resid", "sd")

  # custom list of components.
  if (!is.list(x) || !all(names(x) %in% components))
    stop(paste(c("the conditional probability distribution must be a list with",
      "at least one of the following elements:", components), collapse = " "))
  if (!("coef" %in% names(x)) || !any(c("sd", "resid") %in% names(x)))
    stop("at least the regression coefficients and one between the residulals ",
      "or the residual standard deviation are required.")

  # check the regression coefficients.
  if (!is.null(x$coef)) {

    if (length(x$coef) == 0)
      stop("fitted must a vector of numeric values, ",
        "the regression coefficients for the node given its parents.")
    if (!is.numeric(x$coef) || !all(is.finite(x$coef)))
      stop("fitted must a vector of numeric values, ",
        "the regression coefficients for the node given its parents.")
    if (!is.null(names(x$coef)))
      if (all(names(x$coef) != "(Intercept)"))
        stop("the intercept is missing.")

  }#THEN

  # check the fitted xs.
   if (!is.null(x$fitted)) {

    if (length(x$fitted) == 0)
      stop("coef must a vector of numeric values, ",
        "the fitted values for the node given its parents.")
    if (!is.numeric(x$fitted) || !all(is.finite(x$fitted)))
      stop("coef must a vector of numeric values, ",
        "the fitted values for the node given its parents.")

  }#THEN

 # check the residuals.
  if (!is.null(x$resid)) {

    if (length(x$resid) == 0)
      stop("resid must a vector of numeric values, ",
        "the residuals for the node given its parents.")
    if (!is.numeric(x$resid) || !all(is.finite(x$resid)))
      stop("resid must a vector of numeric values, ",
        "the residuals for the node given its parents.")

  }#THEN

  # check the standard deviation of the residuals.
  if (!is.null(x$sd)) {

    if (!is.non.negative(x$sd))
      stop("sd must be a non-negative number, the standard deviation of the residuals.")

    if (!is.null(x$resid) && !isTRUE(all.equal(x$sd, sd(x$resid))))
      stop("the reported standard deviation of the residuals does not match the observed one.")

  }#THEN

  # one residual for each fitted value.
  if (!is.null(x$resid) && !is.null(x$fitted))
    if (length(x$resid) != length(x$fitted))
      stop("the residuals and the fitted values have different lengths.")

}#CHECK.FIT.GNODE.SPEC

# check one bn.fit.gnode against another.
check.gnode.vs.spec = function(new, old, node) {

  if (is(old, "bn.fit.gnode")) {

    # same number of coefficients.
    if (length(new$coef) != length(old$coefficients))
      stop("wrong number of coefficients for node ", old$node, ".")
    # if the new coefficients have labels, they must match.
    if (!is.null(names(new$coef)))
      check.nodes(names(new$coef), graph = names(old$coefficients),
        min.nodes = length(names(old$coefficients)))
    # same number of residuals.
    if (!is.null(new$resid))
      if (length(new$resid) != length(old$residuals))
        stop("wrong number of residuals for node ", old$node, ".")
    # same number of fitted values.
    if (!is.null(new$fitted))
      if (length(new$fitted) != length(old$fitted.values))
        stop("wrong number of fitted values for node ", old$node, ".")

  }#THEN
  else {

    # add the intercept, which is obviously not among the parents of the node.
    old = c("(Intercept)", old)

    # same number of coefficients.
    if (length(new$coef) != length(old))
      stop("wrong number of coefficients for node ", node, ".")
    # if the new coefficients have labels, they must match.
    if (!is.null(names(new$coef)))
      check.nodes(names(new$coef), graph = old, min.nodes = length(old))

  }#ELSE

}#CHECK.GNODE.VS.SPEC

# check a user-specified bn.fit.dnode.
check.fit.dnode.spec = function(x) {

  # the CPT must be a table.
  if (!is.ndmatrix(x))
    stop("the conditional probability distribution must be a table, ",
      "a matrix or a multidimensional array.")
  # all elements must be probabilities.
  if (!is.probability.vector(x))
    stop("some elements of the conditional probability distributions are ",
      "not probabilities.")

  # convert the CPT into a table object.
  x = as.table(x)

  # compute the dimensions of the table
  dims = dim(x)
  ndims = length(dims)
  # flatten 1xc tables into 1-dimensional ones.
  if ((ndims == 2) && (dims[1] == 1))
    x = flatten.2d.table(x)
  # update dims and ndims.
  dims = dim(x)
  ndims = length(dims)

  # normalization.
  if (ndims == 1) {

    if (sum(x) != 1)
      stop("the probability distribution does not sum to one.")

  }#THEN
  else {

    if (any(margin.table(x,  seq(from = 2, to = ndims)) != 1))
      stop("some conditional probability distributions do not sum to one.")

  }#ELSE

  return(x)

}#CHECK.FIT.DNODE.SPEC

# check one bn.fit.dnode against another.
check.dnode.vs.spec = function(new, old, node, cpt.levels) {

  ndims = length(dim(new))

  if (is(old, "bn.fit.dnode")) {

    # same dimensions.
    if (!identical(dim(new), dim(old$prob)))
      stop("wrong dimensions for node ", old$node, ".")
    # if the new CPT has labels, they must match.
    if (!is.null(dimnames(new)))
      if (!identical(dimnames(new), dimnames(old$prob)))
        stop("wrong levels for node ", old$node, ".")

  }#THEN
  else {

    # same dimensions and labels.
    if (ndims == 1) {

      if (dim(new) != length(cpt.levels[[node]]))
        stop("wrong dimensions for node ", node, ".")
      if (!identical(cpt.levels[[node]], names(new)))
        stop("wrong levels for node ", node, ".")

    }#THEN
    else {

      if (any(dim(new) != sapply(cpt.levels[c(node, old)], length)))
        stop("wrong number of dimensions for node ", node, ".")
      if (!all(names(dimnames(new)) %in% c(node, old)))
        stop("wrong dimensions for node ", node, ".")
      # now that we are sure that the dimensions are the right ones, reorder them
      # to follow the ordering of the parents in the network.
      d = names(dimnames(new))
      new = aperm(new, c(match(node, d), match(old, d)))

      if (!identical(cpt.levels[c(node, old)], dimnames(new)))
        stop("wrong levels for node ", node, ".")

    }#ELSE

  }#ELSE

  return(new)

}#CHECK.DNODE.VS.SPEC

