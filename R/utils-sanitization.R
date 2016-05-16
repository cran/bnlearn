
# is x a vector (as opposed to a matrix)?
is.vector = function(x) {

  is.null(dim(x)) ||
  (length(dim(x)) == 2) &&
  (dim(x)[2] == 1)

}#IS.VECTOR

# is x a real number?
is.real.number = function(x) {

  is.numeric(x) &&
  (length(x) == 1) &&
  is.finite(x)

}#IS.REAL

# is x a vector of real number?
is.real.vector = function(x) {

  is.numeric(x) &&
  all(is.finite(x))

}#IS.REAL.VECTOR

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
is.probability.vector = function(x, zero = FALSE) {

  is.numeric(x) &&
  all(is.finite(x)) &&
  all(x >= 0) &&
  all(x <= 1) &&
  (zero || any(x > 0))

}#IS.PROBABILITY.VECTOR

# is x a single character string?
is.string = function(x) {

  is.character(x) &&
  (length(x) == 1) &&
  !any(is.na(x)) &&
  any(x != "")

}#IS.STRING

# is x a vector of character strings?
is.string.vector = function(x) {

  is.character(x) &&
  !any(is.na(x)) &&
  any(x != "")

}#IS.STRING

is.ndmatrix = function(x) {

  is(x, c("table", "matrix", "array"))

}#IS.NDMATRIX

# is x a symmetric matrix?
is.symmetric = function(x) {

  .Call("is_symmetric",
        matrix = x)

}#IS.SYMMETRIC

# does x satisfy the Cauchy-Schwarz inequality?
is.cauchy.schwarz = function(x) {

  .Call("is_cauchy_schwarz",
        matrix = x)

}#IS.CAUCHY.SCHWARZ

# are all numeric values in the data frame fimite?
check.data.frame.finite = function(x) {

  .Call("data_frame_finite",
        data = x)

}#DATA.FRAME.FINITE

# check the data set.
check.data = function(x, allowed.types = available.data.types,
    allow.levels = FALSE) {

  # check the data are there.
  if (missing(x))
    stop("the data are missing.")
  # x must be a data frame.
  if(!is.data.frame(x))
    stop("the data must be in a data frame.")
  # check the data for NULL/NaN/NA.
  if (missing.data(x))
    stop("the data set contains NULL/NaN/NA values.")
  # check which type of data we are dealing with.
  type = data.type(x)

  # check whether the variables are either all continuous or all discrete.
  if (type %!in% allowed.types)
    stop("valid data types are:\n",
      sprintf("    * %s.\n", data.type.labels[allowed.types]))

  # checks specific to a particular data type.
  if (type %in% discrete.data.types) {

    for (col in names(x)) {

      # check the number of levels of discrete variables, to guarantee that
      # the degrees of freedom of the tests are positive.
      if (nlevels(x[, col]) < 2)
        stop("variable ", col, " must have at least two levels.")

      # warn about levels with zero frequencies, it's not necessarily wrong
      # (data frame subsetting) but sure is fishy.
      if (!allow.levels && any(table(x[, col]) == 0))
        warning("variable ", col, " has levels that are not observed in the data.")

    }#FOR

  }#THEN
  else if (type == "continuous") {

    # all values must be finite to have finite mean and variance.
    check.data.frame.finite(x)

  }#THEN

  return(type)

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
    stop("at most ", max.nodes, " node(s) needed.")
  # minimum number of nodes requirement (usually 1).
  if (length(nodes) < min.nodes)
    stop("at least ", min.nodes, " node(s) needed.")
  # node must be a valid node label.
  if (!is.null(graph)) {

    if (is(graph, "bn")) {

      if (any(nodes %!in% names(graph$nodes)))
        stop("node(s)", paste0(" '", nodes[nodes %!in% names(graph$nodes)], "'"),
             " not present in the graph.")

    }#THEN
    else if (is(graph, "bn.fit")) {

      if (any(nodes %!in% names(graph)))
        stop("node(s)", paste0(" '", nodes[nodes %!in% names(graph)], "'"),
             " not present in the graph.")

    }#THEN
    else if (is.character(graph)) {

      if (any(nodes %!in% graph))
        stop("node(s)", paste0(" '", nodes[nodes %!in% graph], "'"),
             " not present in the graph.")

    }#THEN

  }#THEN

}#CHECK.NODES

# check an arc set.
check.arcs = function(arcs, nodes) {

  # sanitize the set of arcs.
  if (is(arcs, "matrix") || is(arcs, "data.frame")) {

     if (dim(arcs)[2] != 2)
       stop("arc sets must have two columns.")
     if (!all(sapply(arcs, class) == "character"))
       stop("node labels in arc sets must be character strings.")

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
  if (any(arcs %!in% nodes))
    stop("node(s)", paste0(" '", unique(arcs[arcs %!in% nodes]), "'"),
         " not present in the graph.")

  # remove duplicate arcs.
  arcs = unique.arcs(arcs, nodes, warn = TRUE)

  # check there are no loops among the arcs.
  loop = (arcs[, "from"] == arcs[, "to"])

  if (any(loop))
    stop("invalid arcs that are actually loops:\n",
         paste("  ", arcs[loop, 1], "->", arcs[loop, 2], "\n"))

  return(arcs)

}#CHECK.ARCS

# build a valid whitelist.
build.whitelist = function(whitelist, nodes, data, algo, criterion) {

  if (is.null(whitelist)) {

    # no whitelist, nothing to do.
    return(NULL)

  }#THEN

  if (is(whitelist, c("matrix", "data.frame"))) {

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
  if (any(unique(as.vector(whitelist)) %!in% nodes))
    stop("unknown node label present in the whitelist.")
  # check that whitelisted arcs do not violate parametric assumptions.
  whitelist = check.arcs.against.assumptions(whitelist, data, criterion)

  if (algo %in% score.based.algorithms) {

    # the whitelist should contain only directed arcs; extend the implied CPDAG
    # instead of picking arc directions at random to avoid loops.
    whitelist = cpdag.arc.extension(whitelist, nodes = nodes)

  }#THEN
  else if (algo %in% mim.based.algorithms) {

    # all arcs in the whitelist are treated as undirected, because these
    # algorithms operate in the space of undirected graphs.
    whitelist = unique.arcs(arcs.rbind(whitelist, whitelist,
                  reverse2 = TRUE), nodes)

  }#THEN

  # if the whitelist itself contains cycles, no acyclic graph
  # can be learned.
  if (!is.acyclic(whitelist, nodes = nodes,
         directed = (algo %in% c(constraint.based.algorithms, "aracne"))))
    stop("this whitelist does not allow an acyclic graph.")

  return(whitelist)

}#BUILD.WHITELIST

check.arcs.against.assumptions = function(arcs, data, criterion) {

  if (criterion %in% c(available.mixedcg.tests, available.mixedcg.scores)) {

    # arcs cannot point from continuous nodes to discrete nodes.
    if (is.null(arcs)) {

      arcs = list.cg.illegal.arcs(nodes = names(data), variables = data)

    }#THEN
    else {

      arcs = .Call("arcs_cg_assumptions",
                   arcs = arcs,
                   nodes = names(data),
                   data = data)

    }#ELSE

  }#THEN

  return(arcs)

}#CHECK.ARCS.AGAINST.ASSUMPTIONS

list.cg.illegal.arcs = function(nodes, variables) {

  .Call("cg_banned_arcs",
        nodes = nodes,
        variables = variables)

}#LIST.ILLEGAL.ARCS

# build a valid blacklist.
build.blacklist = function(blacklist, whitelist, nodes, algo) {

  if (!is.null(blacklist)) {

    if (is(blacklist, c("matrix", "data.frame"))) {

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
    if (any(unique(as.vector(blacklist)) %!in% nodes))
      stop("unknown node label present in the blacklist.")

    if (algo %in% mim.based.algorithms) {

      # all arcs in the whitelist are treated as undirected, because these
      # algorithms operate in the space of undirected graphs.
      blacklist = arcs.rbind(blacklist, blacklist, reverse2 = TRUE)

    }#THEN

    # drop duplicate rows.
    blacklist = unique.arcs(blacklist, nodes)

  }#THEN

  # update blacklist to agree with whitelist.
  # NOTE: whitelist and blacklist relationship is the same as hosts.allow
  # and hosts.deny.
  if (!is.null(whitelist)) {

    # if x -> y is whitelisted but y -> x is not, it is to be blacklisted.
    to.add = apply(whitelist, 1, function(x)
               is.whitelisted(whitelist, x[c(2, 1)]))
    blacklist = arcs.rbind(blacklist, whitelist[!to.add, c(2, 1)])

    # if x -> y is whitelisted, it is to be removed from the blacklist.
    if (!is.null(blacklist)) {

      blacklist = blacklist[!apply(blacklist, 1,
        function(x){ is.whitelisted(whitelist, x) }),]

      # also drop duplicate rows.
      blacklist = unique.arcs(matrix(blacklist, ncol = 2, byrow = FALSE), nodes)

    }#THEN

  }#THEN

  # set the column names.
  if (!is.null(blacklist))
    colnames(blacklist) = c("from", "to")

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

  # check which type of data we are dealing with.
  type = data.type(data)

  if (!is.null(score)) {

    # check it's a single character string.
    check.string(score)
    # check the score/test label.
    if (score %!in% available.scores)
      stop("valid scores are:\n",
           sprintf("    %-15s %s\n", names(score.labels), score.labels))
    # check if it's the right score for the data (discrete, continuous, mixed).
    if ((type %!in% discrete.data.types) &&
         (score %in% available.discrete.scores))
      stop("score '", score, "' may be used with discrete data only.")
    if ((type != "continuous") && (score %in% available.continuous.scores))
      stop("score '", score, "' may be used with continuous data only.")
    if ((type != "mixed-cg") && (score %in% available.mixedcg.scores))
      stop("score '", score, "' may be used with a mixture of continuous and discrete data only.")

    return(score)

  }#THEN
  else {

    # warn about ordinal data modelled as unordered categorical ones.
    if (type %in% c("ordered", "mixed-do"))
      warning("no score is available for ordinal data, disregarding the ordering of the levels.")

    if (type %in% discrete.data.types)
      return("bic")
    else if (type == "continuous")
      return("bic-g")
    else if (type == "mixed-cg")
      return("bic-cg")

  }#ELSE

}#CHECK.SCORE

# check whether a score is score equivalent.
is.score.equivalent = function(score, nodes, extra) {

  # log-likelihood for discrete and Gaussian data is always score equivalent.
  if (score %in% c("loglik", "loglik-g"))
    return(TRUE)
  # same with AIC and BIC.
  if (score %in% c("aic", "aic-g", "bic", "bic-g"))
    return(TRUE)
  # BDe and BGe are score equivalent if they have uniform priors (e.g. BDeu and BGeu).
  else if ((score %in% c("bde", "bge")) && (extra$prior == "uniform"))
    return(TRUE)

  # a conservative default.
  return(FALSE)

}#IS.SCORE.EQUIVALENT

# check whether a score is decomposable.
is.score.decomposable = function(score, nodes, extra) {

  # Castelo & Siebes prior is not decomposable.
  if ((score %in% c("bde", "bds", "mbde", "bge")) && (extra$prior == "cs"))
    return(FALSE)

  # a sensible default.
  return(TRUE)

}#IS.SCORE.DECOMPOSABLE

# check test labels.
check.test = function(test, data) {

  # check which type of data we are dealing with.
  type = data.type(data)

  if (!missing(test) && !is.null(test)) {

    # check it's a single character string.
    check.string(test)
    # check the score/test label.
    orig.length = unlist(options("warning.length" = 4000))
    if (test %!in% available.tests)
      stop("valid tests are:\n", sprintf("    %-15s %s\n",
           names(test.labels), test.labels))
    # check if it's the right test for the data (discrete, continuous).
    if ((type != "ordered") && (test %in% available.ordinal.tests))
      stop("test '", test, "' may be used with ordinal data only.")
    if ((type %!in% discrete.data.types) && (test %in% available.discrete.tests))
      stop("test '", test, "' may be used with discrete data only.")
    if ((type != "continuous") && (test %in% available.continuous.tests))
      stop("test '", test, "' may be used with continuous data only.")
    if ((type != "mixed-cg") && (test %in% available.mixedcg.tests))
      stop("test '", test, "' may be used with a mixture of continuous and discrete data only.")
    options("warning.length" = orig.length)

    return(test)

  }#THEN
  else {

    if (type == "ordered")
      return("jt")
    else if (type %in% c("factor", "mixed-do"))
      return("mi")
    else if (type == "continuous")
      return("cor")
    else if (type == "mixed-cg")
      return("mi-cg")

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
    stop("valid tests are:\n",
         sprintf("    %-15s %s\n", names(test.labels), test.labels),
         "  valid scores are:\n",
         sprintf("    %-15s %s\n", names(score.labels), score.labels))

  return(criterion)

}#CHECK.CRITERION

# check loss functions' labels.
check.loss = function(loss, data, bn) {

  # check which type of data we are dealing with.
  type = data.type(data)

  if (!is.null(loss)) {

    # check it's a single character string.
    check.string(loss)
    # check the score/test label.
    if (loss %!in% loss.functions)
      stop("valid loss functions are:\n",
           sprintf("    %-15s %s\n", names(loss.labels), loss.labels))
    if ((type %!in% discrete.data.types) && (loss %in% discrete.loss.functions))
      stop("loss function '", loss, "' may be used with discrete data only.")
    if ((type != "continuous") && (loss %in% continuous.loss.functions))
      stop("loss function '", loss, "' may be used with continuous data only.")
    if ((type != "mixed-cg") && (loss %in% mixedcg.loss.functions))
      stop("loss function '", loss, "' may be used with a mixture of continuous and discrete data only.")

    return(loss)

  }#THEN
  else {

    if ((is.character(bn) && (bn %in% classifiers)) ||
         is(bn, c("bn.naive", "bn.tan")))
      return("pred")
    if (type %in% discrete.data.types)
      return("logl")
    else if (type == "continuous")
      return("logl-g")
    else if (type == "mixed-cg")
      return("logl-cg")

  }#ELSE

}#CHECK.LOSS

# check the method used to fit the parameters of the network.
check.fitting.method = function(method, data) {

  # check which type of data we are dealing with.
  type = data.type(data)

  if (!is.null(method)) {

    # check it's a single character string.
    check.string(method)
    # check the score/test label.
    if (method %!in% available.fitting.methods)
      stop("valid fitting methods are:\n",
           sprintf("    %-15s %s\n", names(fitting.labels), fitting.labels))
    # bayesian parameter estimation is implemented only for discrete data.
    if ((type %in% c("continuous", "mixed-cg")) && (method == "bayes"))
      stop("Bayesian parameter estimation for (conditional) Gaussian Bayesian networks is not implemented.")

    return(method)

  }#THEN
  else {

    return("mle")

  }#ELSE

}#CHECK.FITTING.METHOD

# check the method used for prediction.
check.prediction.method = function(method, data) {

  if (!is.null(method)) {

    # check it's a single character string.
    check.string(method)
    # check the score/test label.
    if (method %!in% available.prediction.methods)
      stop("valid prediction methods are:\n",
        sprintf("    %-15s %s\n", names(prediction.labels), prediction.labels))

    return(method)

  }#THEN
  else {

    return("parents")

  }#ELSE

}#CHECK.PREDICTION.METHOD

# check the method used to discretize the data.
check.discretization.method = function(method) {

  if (!is.null(method)) {

    # check it's a single character string.
    check.string(method)
    # check the score/test label.
    if (method %!in% available.discretization.methods)
      stop("valid discretization methods are:\n",
           sprintf("    %-15s %s\n", names(discretization.labels),
                   discretization.labels))

    return(method)

  }#THEN
  else {

      return("quantile")

  }#ELSE

}#CHECK.DISCRETIZATION.METHOD

# check the estimator for the mutual information.
check.mi.estimator = function(estimator, data) {

  # check which type of data we are dealing with.
  type = data.type(data)

  if (!is.null(estimator)) {

    # check it's a single character string.
    check.string(estimator)
    # check the score/estimator label.
    if (estimator %!in% available.mi)
      stop("valid estimators are:\n",
           sprintf("    %-15s %s\n", names(mi.estimator.labels),
                   mi.estimator.labels))
    # check if it's the right estimator for the data (discrete, continuous).
    if ((type %!in% discrete.data.types) &&
        (estimator %in% available.discrete.mi))
      stop("estimator '", estimator, "' may be used with discrete data only.")
    if ((type != "continuous") && (estimator %in% available.continuous.mi))
      stop("estimator '", estimator, "' may be used with continuous data only.")

    return(estimator)

  }#THEN
  else {

    if (type %in% discrete.data.types)
      return("mi")
    else
      return("mi-g")

  }#ELSE

}#CHECK.MI.ESTIMATOR

# is the data of a particular type?
data.type = function(data) {

  .Call("data_type",
        data = data)

}#DATA.TYPE

# there are missing data?
missing.data = function(data) {

  !all(complete.cases(data))

}#MISSING.DATA

# check the imaginary sample size.
check.iss = function(iss, network, data) {

  # check which type of data we are dealing with.
  type = data.type(data)

  if (!is.null(iss)) {

    # validate the imaginary sample size.
    if (!is.positive(iss) || (iss < 1))
      stop("the imaginary sample size must be a numeric value greater than 1.")
    # if iss = 1 the bge is NaN, if iss = 2 and phi = "heckerman" the
    # computation stops with the following error:
    # Error in solve.default(phi[A, A]) :
    #   Lapack routine dgesv: system is exactly singular
    if((type == "continuous") && (iss < 3))
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

    if (phi %!in% c("heckerman", "bottcher"))
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
    if (any(names(exp) %!in% names(data)) || (length(names(exp)) == 0))
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

    # check whether there is a penalization coefficient stored in the bn object,
    # use the default for the score function otherwise.
    if (!is.null(network$learning$args$k) && (score == network$learning$test))
      k = network$learning$args$k
    else
      k = ifelse((score %in% c("aic", "aic-g")), 1, log(nrow(data))/2)

  }#ELSE

  return(k)

}#CHECK.PENALTY

# sanitize prior distributions over the graph space.
check.graph.prior = function(prior, network) {

  if (is.null(prior)) {

    # check whether there is a graph prior stored in the bn object, use the
    # uniform one otherwise.
    if (!is.null(network$learning$args$prior))
      prior = network$learning$args$prior
    else
      prior = "uniform"

  }#THEN
  else {

    # check whether prior is a string.
    check.string(prior)
    # check whether the label matches a known prior.
    if (prior %!in% prior.distributions)
      stop("valid prior distributions are: ",
        paste(prior.distributions, collapse = " "), ".")

  }#ELSE

  return(prior)

}#CHECK.GRAPH.PRIOR

# check the sparsity parameter of the prior distribution over the graph space.
check.graph.sparsity = function(beta, prior, network, data, learning = FALSE) {

  default.beta =
    list("uniform" = NULL, "vsp" = 1/ncol(data),
      "cs" = cs.completed.prior(data.frame(character(0), character(0),
             numeric(0)), names(data)))

  if (is.null(beta)) {

    # check whether there is a graph prior stored in the bn object, use the
    # uniform one otherwise.
    if (!is.null(network$learning$args$prior))
      beta = network$learning$args$beta
    else
      beta = default.beta[[prior]]

  }#THEN
  else {

    if (prior == "uniform") {

      warning("unused argument beta.")
      beta = NULL

    }#THEN
    else if (prior == "vsp") {

      if (!is.probability(beta) || (beta >= 1))
       stop("beta must be a probability smaller than 1.")

    }#THEN
    else if (prior == "cs") {

      # arcs' prior probabilities should be provided in a data frame.
      if (!is.data.frame(beta) || (ncol(beta) != 3) ||
          !identical(colnames(beta), c("from", "to", "prob")))
        stop("beta must be a data frame with three colums: 'from', 'to' and 'prob'.")
      # the probs cloumns must contain only probabilities.
      if (!is.probability.vector(beta$prob, zero = TRUE))
        stop("arcs prior must contain only probabilities.")
      # check that the first two columns contain only valid arcs.
      check.arcs(beta[, c("from", "to")], nodes = names(data))

      # complete the user-specified prior.
      beta = cs.completed.prior(beta, names(data), learning)

    }#THEN

  }#ELSE

  return(beta)

}#CHECK.GRAPH.SPARSITY

check.maxp = function(maxp, data) {

  if (is.null(maxp)) {

    maxp = Inf

  }#THEN
  else if (!isTRUE(all.equal(maxp, Inf))) {

    if (!is.positive.integer(maxp))
      stop("maxp must be a positive number.")
    if (maxp >= ncol(data))
      warning("maximum number of parents should be lower than the number of nodes, the limit will be ignored.")

  }#ELSE

  return(as.numeric(maxp))

}#CHECK.MAXP

# sanitize the extra arguments passed to the network scores.
check.score.args = function(score, network, data, extra.args, learning = FALSE) {

  # check the imaginary sample size.
  if (score %in% c("bde", "bds", "mbde", "bge"))
    extra.args$iss = check.iss(iss = extra.args$iss,
      network = network, data = data)

  # check the graph prior distribution.
  if (score %in% c("bde", "bds", "mbde", "bge"))
    extra.args$prior = check.graph.prior(prior = extra.args$prior,
      network = network)

  # check the sparsity parameter of the graph prior distribution.
  if (score %in% c("bde", "bds", "mbde", "bge"))
    extra.args$beta = check.graph.sparsity(beta = extra.args$beta,
      prior = extra.args$prior, network = network, data = data,
      learning = learning)

  # check the list of the experimental observations in the data set.
  if (score == "mbde")
    extra.args$exp = check.experimental(exp = extra.args$exp,
      network = network, data = data)

  # check the likelihood penalty.
  if (score %in% c("aic", "bic", "aic-g", "bic-g", "aic-cg", "bic-cg"))
    extra.args$k = check.penalty(k = extra.args$k, network = network,
      data = data, score = score)

  # check phi estimator.
  if (score == "bge")
    extra.args$phi = check.phi(phi = extra.args$phi,
      network = network, data = data)

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
check.cpq.args = function(fitted, extra.args, method, action) {

  if (method %in% c("ls", "lw")) {

    if (!is.null(extra.args$n)) {

      if (!is.positive.integer(extra.args$n))
        stop("the number of observations to be sampled must be a positive integer number.")

    }#THEN
    else {

      # this is a rule of thumb, the error of the estimate has no closed-form
      # expression (Koller & Friedman).
      if (!is(fitted, "bn.fit.gnet"))
        extra.args$n = 5000 * max(1, round(log10(nparams.fitted(fitted))))
      else
        extra.args$n = 500 * nparams.fitted(fitted)

    }#ELSE

    if (!is.null(extra.args$batch)) {

      if ((action == "cpdist") && (method == "lw")) {

        extra.args$batch = NULL
        warning(" 'batch' will be ignored for speed and memory efficience.")

      }#THEN
      else {

        if (!is.positive.integer(extra.args$batch))
          stop("the number of observations to be sampled must be a positive integer number.")

        if (extra.args$batch > extra.args$n) {

          warning("cannot generate a batch bigger than the whole generated data set.")
          warning("batch size will be ignored.")

        }#THEN

      }#ELSE

    }#THEN
    else {

      # perform small simulations in a single batch, and split larger ones.
      extra.args$batch = min(extra.args$n, 10^4)

    }#ELSE

    if (!is.null(extra.args$query.nodes)) {

      check.nodes(extra.args$query.nodes, graph = fitted)

    }#THEN

  }#THEN

  # warn about unused arguments.
  check.unused.args(extra.args, cpq.extra.args[[method]])

  return(extra.args)

}#CHECK.CPQ.ARGS

# sanitize the extra arguments passed to loss functions.
check.loss.args = function(loss, bn, nodes, data, extra.args) {

  valid.args = loss.extra.args[[loss]]

  if (loss %in% c("pred", "pred-lw", "pred-lw-cg", "cor", "cor-lw", "cor-lw-cg",
                  "mse", "mse-lw", "mse-lw-cg")) {

    if (!is.null(extra.args$target)) {

      if (!is.string(extra.args$target) || (extra.args$target %!in% nodes))
        stop("target node must be a single, valid node label for the network.")

      # in hybrid networks, check the target has the right data type.
      if (loss %in% c("cor-lw-cg", "mse-lw-cg"))
        if (!is(data[, extra.args$target], "numeric"))
          stop("the target node must be a continuous variable.")
      if (loss == "pred-lw-cg")
        if (!is(data[, extra.args$target], "factor"))
          stop("the target node must be a factor.")

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
    if ((bn %in% classifiers) || is(bn, c("bn.naive", "bn.tan"))) {

      extra.args$prior = check.classifier.prior(extra.args$prior, data[, extra.args$target])
      valid.args = c(valid.args, "prior")

    }#THEN

  }#THEN

  if (loss %in% c("pred-lw", "pred-lw-cg", "cor-lw", "cor-lw-cg", "mse-lw",
                  "mse-lw-cg")) {

    # number of particles for likelihood weighting.
    if (!is.null(extra.args$n)) {

      if (!is.positive.integer(extra.args$n))
        stop("the number of observations to be sampled must be a positive integer number.")

    }#THEN
    else {

      extra.args$n = 500

    }#ELSE

    # which nodes to predict from.
    if (!is.null(extra.args$from))
      check.nodes(extra.args$from, graph = names(data), min.nodes = 1)
    else
      extra.args$from = setdiff(names(data), extra.args$target)

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

  # check which type of data we are dealing with.
  type = data.type(data)

  if (method == "hartemink") {

    if (type %in% discrete.data.types) {

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
        if (idisc %!in% other.methods)
          stop("valid initial discretization methods are:\n",
               sprintf("    %-15s %s\n", other.methods,
                       discretization.labels[other.methods]))

      }#THEN
      else {

        # default to quantile discretization as per Hartemink's recommendation.
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

  if (method == "tree.bayes") {

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

}#CHECK.CLASSIFIER.ARGS

# warn about unused arguments.
check.unused.args = function(dots, used.args) {

  unused.args = (names(dots) %!in% used.args)
  if (any(unused.args))
    warning("unused argument(s):", paste0(" '", names(dots)[unused.args], "'"), ".")

}#CHECK.UNUSED.ARGS

# take care of meaningless dots arguments in plot functions.
sanitize.plot.dots = function(dots, meaningless) {

  # warn about them.
  if (any(names(dots) %in% meaningless))
    warning("arguments ", paste(meaningless, collapse = ", "),
      " will be silently ignored.")
  # nuke them from orbit.
  for (m in meaningless)
    dots[[m]] = NULL

  return(dots)

}#PROCESS.PLOT.DOTS

# check the the target nominal type I error rate
check.alpha = function(alpha, network = NULL) {

  # check the the target nominal type I error rate
  if (!missing(alpha) && !is.null(alpha)) {

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

    if (!missing(B) && !is.null(B)) {

      if (!is.positive.integer(B))
        stop("the number of permutations/bootstrap replications must be a positive integer number.")

      B = as.integer(B)

    }#THEN
    else {

      if (criterion %in% semiparametric.tests)
        B = 100L
      else
        B = 5000L

    }#ELSE

  }#THEN
  else {

    if (!missing(B) && !is.null(B))
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
    if (any(colnames(amat) %!in% nodes))
      stop("node (column label) not present in the graph.")
  # column names must be valid node labels.
  if (!is.null(rownames(amat)))
    if (any(rownames(amat) %!in% nodes))
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

# check an object of class bn.fit.
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
  if (ncol(strength) %!in% 3:4)
    stop("objects of class 'bn.strength' must have 3 or 4 columns.")
  if (!identical(names(strength), c("from", "to", "strength")) &&
      !identical(names(strength), c("from", "to", "strength", "direction")))
    stop("objects of class 'bn.strength' must be data frames with column names ",
         "'from', 'to', 'strength' and (optionally) 'direction'.")
  if (any(c("method", "threshold") %!in% names(attributes(strength))))
    stop("objects of class 'bn.strength' must have a 'method' and a 'strength' attributes.")
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

  # check which type of data we are dealing with.
  type = data.type(data)

  # the number of variables must be the same
  if (length(names(bn$nodes)) != ncol(data))
    stop("the network and the data have different numbers of variables.")
  # the variables must be the same.
  if (length(setdiff(names(bn$nodes), names(data))) != 0)
    stop("the variables in the data and in the network do not match.")
  # data type versus network structure.
  if (type == "mixed-cg")
    check.arcs.against.assumptions(bn$arcs, data, "mi-cg")

}#CHECK.BN.VS.DATA

# check bn.fit metadata against the data it's used with.
check.fit.vs.data = function(fitted, data, subset) {

  fitted.names = names(fitted)
  # check which type of data we are dealing with.
  dtype = data.type(data)

  if (missing(subset)) {

    # the number of variables must be the same.
    if (length(fitted.names) != ncol(data))
      stop("the network and the data have different numbers of variables.")
    # the variables must be the same.
    if (length(setdiff(fitted.names , names(data))) != 0)
      stop("the variables in the data and in the network do not match.")

    subset = fitted.names

  }#THEN
  else {

    # the number of variables must not exceed that of the network.
    if (length(subset) > length(fitted.names))
      stop("the data have more variables than the network.")
    # all the variables in the subset must be present in the data.
    absent = (subset %!in% names(data))
    if (any(absent))
      stop("required variables '", paste(subset[absent], collapse = " "),
           "' are not present in the data.")

  }#ELSE

  .Call("fitted_vs_data",
        fitted = fitted,
        data = data,
        subset = subset)

}#CHECK.FIT.VS.DATA

# check bn.fit.{d,g}node metadata against the data it's used with.
check.fit.node.vs.data = function(fitted, data) {

  relevant = c(fitted$node, fitted$parents)
  # check which type of data we are dealing with.
  type = data.type(data)

  # check whether all relevant nodes are in the data.
  if (any(relevant %!in% names(data)))
    stop("not all required nodes are present in the data.")
  # data type versus network type.
  if (is(fitted, "bn.fit.dnode") && (type == "continuous"))
      stop("continuous data and discrete network.")
  if (is(fitted, "bn.fit.gnode") &&
      (type %in% discrete.data.types))
    stop("discrete data and continuous network.")
  # double-check the levels of the variables against those of the nodes.
  if (is(fitted, "bn.fit.dnode")) {

    for (node in relevant) {

      data.levels = levels(data[, node])
      if (length(relevant) == 1)
        node.levels = dimnames(fitted$prob)[[1]]
      else
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

  if ((lty %!in% 0:6) && (lty %!in% lty.strings))
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

  if (algorithm %!in% ok)
       stop("valid learning algorithms are:\n",
            sprintf("    %-15s %s\n", ok, method.labels[ok]))

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
      if ("test" %!in% names(args))
        if (bn$learning$test %in% available.tests)
          args$test = bn$learning$test

      # set the appropriate value for the optimization flag.
      if ("optimized" %!in% names(args))
        args$optimized = bn$learning$optimized

      # pass along all the parameters in bn$learning$args.
      if (length(bn$learning$args) > 0) {

        if ("alpha" %!in% names(args))
          args$alpha = bn$learning$args$alpha

        if ("test" %in% names(args) && ("B" %!in% names(args)))
          if (args$test %in% resampling.tests)
            args$B = bn$learning$args$B

      }#THEN

    }#THEN
    else if (algorithm %in% score.based.algorithms) {

      if ("score" %!in% names(args))
        if (bn$learning$test %in% available.scores)
          args$score = bn$learning$test

      # set the appropriate value for the optimization flag.
      if ("optimized" %!in% names(args))
        args$optimized = bn$learning$optimized

      # pass along the relevant parameters in bn$learning$args if the score
      # function is the same (hint: different scores have paramenters with
      # the same name but different meanings).
      if (("score" %in% names(args)) && (args$score == bn$learning$test))
        for (arg in names(bn$learning$args))
          if ((arg %!in% names(args)) && (arg %in% (score.extra.args[[args$score]])))
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

# check a prior distribution against the observed variable.
check.classifier.prior = function(prior, training) {

  if (missing(prior) || is.null(prior)) {

    # use the empirical probabilities in the fitted network, or a flat prior
    # as a last resort.
    if (is(training, c("bn.fit.dnode", "bn.fit.onode")))
      prior = training$prob
    else
      prior = rep(1, nlevels(training))

  }#THEN
  else {

    if (is(training, c("bn.fit.dnode", "bn.fit.onode")))
      nlvls = dim(training$prob)[1]
    else
      nlvls = nlevels(training)

    if (length(prior) != nlvls)
      stop("the prior distribution and the training variable have a different number of levels.")
    if (!is.nonnegative.vector(prior))
      stop("the prior distribution must be expressed as a probability vector.")

    # make sure the prior probabilities sum to one.
    prior = prior / sum(prior)

  }#ELSE

  return(prior)

}#CHECK.CLASSIFIER.PRIOR

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
check.fit.gnode.spec = function(x, node) {

  components =  c("coef", "fitted", "resid", "sd", "configs")
  labels = c(fitted = "fitted values", resid = "residuals")

  # custom list of components.
  if (!is.list(x) || any(names(x) %!in% components))
    stop("the conditional probability distribution for node ", node,
         " must be a list with at least one of the following elements:",
         paste0(" '", components, "'"), ".")
  if (("coef" %!in% names(x)) || !any(c("sd", "resid") %in% names(x)))
    stop("at least the regression coefficients and either the residuals or ",
      "the residual standard deviation are required for node ", node, ".")
  # discrete parents' configurations only belong to bn.fit.cgnode.
  if ((length(dim(x$coef)) != 2) && ("configs" %in% names(x)))
    stop("no parents' configurations are needed with a single set of coefficients.")

  if (!is.null(x$coef))
    if ((length(x$coef) == 0) || !is.real.vector(x$coef))
      stop("coef must be a vector or a matrix of numeric values, the ",
        "regression coefficients for node ", node, " given its parents.")

  for (comp in c("fitted", "resid"))
    if (!is.null(x[[comp]]))
      if ((length(x[[comp]]) == 0) || !is.real.vector(x[[comp]]))
        stop(comp, " must be a vector of numeric values, the ",
          labels[comp], " for node ", node, " given its parents.")

  if (!is.null(x$configs))
    if ((length(x$configs) == 0) || !is.factor(x$configs))
      stop("the discrete parents' configurations for node ", node,
        " must be a factor.")

  # check the standard deviation of the residuals.
  if (!is.null(x$sd)) {

    if (!is.nonnegative.vector(x$sd))
      stop("sd must be a non-negative number, the standard deviation of the ",
        "residuals of node ", node, ".")

    if ((length(dim(x$coef)) == 2) && (length(x$sd) != ncol(x$coef)) ||
        (length(dim(x$coef)) == 1) && (length(x$sd) > 1) ||
        (length(dim(x$sd)) > 1))
      stop("the dimensions of sd and coef do not match.")

    if (!is.null(x$resid) && (length(x$sd) == 1)) {

      adj.sd = cgsd(x$resid, p = length(x$coef))

      if (!isTRUE(all.equal(x$sd, adj.sd, check.attributes = FALSE, tol = 0.0005)))
        stop("the reported standard deviation of the residuals of node ", node,
          " does not match the observed one.")

    }#THEN
    else if (!is.null(x$resid) && (length(x$sd) > 1) && !is.null(x$configs)) {

      adj.sd = cgsd(x$resid, configs = x$configs, p = nrow(x$coef))

      if (!isTRUE(all.equal(x$sd, adj.sd, check.attributes = FALSE, tol = 0.0005)))
        stop("the reported standard deviation of the residuals of node ", node,
          " does not match the observed one.")

    }#THEN

  }#THEN
  else if (!is.null(x$resid)) {

    # compute the standard error from the residuals if possible.
    if ((length(dim(x$coef)) == 2) && is.null(x$config))
      stop("sd is missing, and parents' configurations are required to compute it.")

    x$sd = cgsd(x$resid, configs = x$configs,
             p = ifelse(is.matrix(x$coef), nrow(x$coef), length(x$coef)))

  }#THEN

  # one residual for each fitted value.
  if (!is.null(x$resid) && !is.null(x$fitted))
    if (length(x$resid) != length(x$fitted))
      stop("the residuals and the fitted values of node ", node,
        " have different lengths.")
  # if any, one parent configuration for each fitted value and residual.
  if (!is.null(x$config)) {

    if (!is.null(x$resid) && (length(x$configs) != length(x$resid)))
      stop("parents' configurations and residuals of node ", node,
        " have different lengths.")
    if (!is.null(x$fitted) && (length(x$configs) != length(x$fitted)))
      stop("parents' configurations and fitted values of node ", node,
        " have different lengths.")

  }#THEN

  return(x)

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
    else
      names(new$coef) = names(old$coef)
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
    else
      names(new$coef) = old

  }#ELSE

  return(new)

}#CHECK.GNODE.VS.SPEC

# check one bn.fit.gnode against another.
check.cgnode.vs.spec = function(new, old, node) {

  if (is(old, "bn.fit.cgnode")) {

    # right dimensions for the coefficients.
    if (!is(new$coef, "matrix") ||
        !identical(dim(new$coef), dim(old$coefficients)))
      stop("the regression coefficients for node ", old$node, " must be ",
           "in a matrix with one column for each discrete parent and one ",
           "coefficient for each continuous parent.")
    # if the new coefficients have labels, they must match.
    if (!is.null(rownames(new$coef)))
      check.nodes(rownames(new$coef), graph = rownames(old$coefficients),
        min.nodes = length(rownames(old$coefficients)))
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

    # nothing useful can possibly be checked from the parent labels.

  }#ELSE

}#CHECK.CGNODE.SPEC

# check a user-specified bn.fit.dnode.
check.fit.dnode.spec = function(x, node) {

  # the CPT must be a table.
  if (!is.ndmatrix(x))
    stop("the conditional probability distribution of node ", node,
      " must be a table, a matrix or a multidimensional array.")
  # all elements must be probabilities.
  if (!is.probability.vector(x))
    stop("some elements of the conditional probability distributions of node ",
      node, " are not probabilities.")

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

    if (abs(sum(x) - 1) > 0.01)
      stop("the probability distribution of node ", node, " does not sum to one.")

    x = x / sum(x)

  }#THEN
  else {

    check.sum = sapply(margin.table(x,  seq(from = 2, to = ndims)),
                  function(x) abs(sum(x) - 1) > 0.01)

    if (any(check.sum))
      stop("some conditional probability distributions of node ", node, " do not sum to one.")

    x = prop.table(x, margin = seq(ndims)[-1])

  }#ELSE

  return(x)

}#CHECK.FIT.DNODE.SPEC

# check one bn.fit.dnode against another.
check.dnode.vs.spec = function(new, old, node, cpt.levels) {

  ndims = length(dim(new))

  if (is(old, c("bn.fit.dnode", "bn.fit.onode"))) {

    # same dimensions.
    if (!identical(dim(new), dim(old$prob)))
      stop("wrong dimensions for node ", old$node, ".")
    # if the new CPT has labels, they must match (and be in the same order too).
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

      # check whether all the CPT inclues all the relevant variables.
      if(!setequal(names(dimnames(new)), c(node, old)))
        stop("wrong dimensions for node ", node, ".")
      # now that we are sure that the dimensions are the right ones, reorder
      # them to follow the ordering of the parents in the network.
      d = names(dimnames(new))
      new = aperm(new, c(match(node, d), match(old, d)))
      # check whether the margins of the CPT are right.
      if (any(dim(new) != sapply(cpt.levels[c(node, old)], length)))
        stop("wrong number of dimensions for node ", node, ".")
      if (!identical(cpt.levels[c(node, old)], dimnames(new)))
        stop("wrong levels for node ", node, ".")

    }#ELSE

  }#ELSE

  return(new)

}#CHECK.DNODE.VS.SPEC

# check evidence in list format for mutilated networks.
check.mutilated.evidence = function(evidence, graph) {

  # check whether evidence is there.
  if (missing(evidence))
    stop("evidence must be a list with elements named after the nodes in the graph.")
  # if evindence is TRUE there's nothing to check.
  if (identical(evidence, TRUE))
    return(TRUE)
  # check whether evidence is a named list.
  if(!is(evidence, "list"))
    stop("evidence must be a list with elements named after the nodes in the graph.")
  # check the node labels in evidence.
  check.nodes(names(evidence), graph = graph)

  if (is(graph, "bn")) {

    # check the network is completely directed.
    if (!is.dag(graph$arcs, names(graph$nodes)))
      stop("the graph is only partially directed.")

  }#THEN
  else {

     # check the evidence is appropriate for the nodes.
     for (fixed in names(evidence)) {

       # extract the node and the evidence.
       cur = graph[[fixed]]
       ev = evidence[[fixed]]

       if (is(cur, c("bn.fit.dnode", "bn.fit.onode"))) {

         if (is.factor(ev))
           evidence[[fixed]] = ev = as.character(ev)

         if (!is.string.vector(ev) || any(ev %!in% dimnames(cur$prob)[[1]]))
           stop("the evidence for node ", fixed, " must be valid levels.")

       }#THEN
       else if (is(cur, "bn.fit.gnode")) {

         # for continuous nodes evidence must be real numbers.
         if (!is.real.vector(ev) || (length(ev) %!in% 1:2))
           stop("the evidence ", fixed, " must be a real number or a finite interval.")
         storage.mode(ev) = "double"
         # amke sure interval boundaries are in the right order.
         evidence[[fixed]] = sort(ev)

       }#THEN

     }#FOR

  }#THEN

  return(evidence)

}#CHECK.MUTILATED.EVIDENCE

