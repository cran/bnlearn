
# check score labels.
check.score = function(score, data, allowed = available.scores) {

  # check which type of data we are dealing with.
  type = data.type(data)

  if (!is.null(score)) {

    # check the score label.
    check.label(score, choices = allowed, labels = score.labels,
      argname = "score", see = "bnlearn-package")
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
  if (score %in% c("loglik", "loglik-g", "pred-loglik", "pred-loglik-g"))
    return(TRUE)
  # same with AIC and BIC.
  if (score %in% c("aic", "aic-g", "bic", "bic-g"))
    return(TRUE)
  # BDe and BGe can score equivalent depending on the graph prior.
  else if ((score %in% c("bde", "bge")) &&
           (extra$prior %in% c("uniform", "marginal", "vsp")))
    return(TRUE)

  # a conservative default (BDla, BDs, BDj, all *-cg scores).
  return(FALSE)

}#IS.SCORE.EQUIVALENT

# check whether a score is decomposable.
is.score.decomposable = function(score, extra) {

  # Castelo & Siebes prior is not decomposable.
  if ((score %in% c("bde", "bds", "bdj", "mbde", "bdla", "bge")) &&
      (extra$prior %in% c("cs", "marginal")))
    return(FALSE)

  # a sensible default.
  return(TRUE)

}#IS.SCORE.DECOMPOSABLE

# sanitize the extra arguments passed to the network scores.
check.score.args = function(score, network, data, extra.args, learning = FALSE) {

  # check the imaginary sample size.
  if (has.argument(score, "iss", score.extra.args))
    extra.args$iss = check.iss(iss = extra.args$iss,
      network = network)

  # check the graph prior distribution.
  if (has.argument(score, "prior", score.extra.args))
    extra.args$prior = check.graph.prior(prior = extra.args$prior,
      network = network)

  # check the sparsity parameter of the graph prior distribution.
  if (has.argument(score, "beta", score.extra.args))
    extra.args$beta = check.graph.hyperparameters(beta = extra.args$beta,
      prior = extra.args$prior, network = network, data = data,
      learning = learning)

  # check the list of the experimental observations in the data set.
  if (has.argument(score, "exp", score.extra.args))
    extra.args$exp = check.experimental(exp = extra.args$exp,
      network = network, data = data)

  # check the likelihood penalty.
  if (has.argument(score, "k", score.extra.args))
    extra.args$k = check.penalty(k = extra.args$k, network = network,
      data = data, score = score)

  # check the number of scores to average.
  if (has.argument(score, "l", score.extra.args))
    extra.args$l = check.l(l = extra.args$l)

  # check the normal-Wishart prior arguments.
  if (has.argument(score, "nu", score.extra.args))
    extra.args$nu = check.nu(nu = extra.args$nu, network = network, data = data)

  if (has.argument(score, "iss.mu", score.extra.args))
    extra.args$iss.mu = check.iss.mu(iss.mu = extra.args$iss.mu,
      network = network)

  if (has.argument(score, "iss.w", score.extra.args))
    extra.args$iss.w = check.iss.w(iss.w = extra.args$iss.w,
      network = network, data = data)

  # check the test data for predictive scores.
  if (has.argument(score, "newdata", score.extra.args))
    extra.args$newdata = check.newdata(newdata = extra.args$newdata,
      network = network, data = data)

  # check the R function implementing a custom score.
  if (has.argument(score, "fun", score.extra.args))
    extra.args$fun = check.custom.score.function(fun = extra.args$fun)
  if (has.argument(score, "args", score.extra.args))
    extra.args$args = check.custom.score.arguments(args = extra.args$args)

  check.unused.args(extra.args, score.extra.args[[score]])

  return(extra.args)

}#CHECK.SCORE.ARGS

# check the imaginary sample size of the Dirichlet prior.
check.iss = function(iss, network) {

  if (!is.null(iss)) {

    # validate the imaginary sample size.
    if (!is.positive(iss))
      stop("the imaginary sample size must be a positive number.")
    if (iss < 1)
      warning("very small imaginary sample size, the results may be affected by numeric problems.")

  }#THEN
  else {

    # check whether there is an imaginary sample size stored in the bn object;
    # otherwise use a the de facto standard value of 1.
    if (!is.null(network$learning$args$iss))
      iss = network$learning$args$iss
    else
      iss = 1

  }#ELSE

  # coerce iss to integer.
  return(as.numeric(iss))

}#CHECK.ISS

# check the prior mean in the BGe score.
check.nu = function(nu, network, data) {

  if (!is.null(nu)) {

    # check type and length.
    if (!is.real.vector(nu))
      stop("the prior mean vector must contain real numbers.")
    if (length(nu) != ncol(data))
      stop("the length of the prior mean vector does not match the number of variables in the data.")

    # check the names.
    if (is.null(names(nu)))
      names(nu) = names(data)
    else if (!setequal(names(nu), names(data)))
      stop("the prior mean vector and the data contain variables with different names.")

  }#THEN
  else {

    if (!is.null(network$learning$args$nu))
      nu = network$learning$args$nu
    else
      nu = structure(colMeans(data), names = names(data))

  }#ELSE

  return(nu)

}#CHECK.NU

# check the imaginary sample size for the normal prior over the mean.
check.iss.mu = function(iss.mu, network) {

  if (!is.null(iss.mu)) {

    # validate the imaginary sample size.
    if (!is.positive(iss.mu))
      stop("the imaginary sample size must be a positive number.")
    if (iss.mu < 1)
      warning("very small imaginary sample size, the results may be affected by numeric problems.")

  }#THEN
  else {

    # check whether there is an imaginary sample size stored in the bn object;
    # otherwise use a the de facto standard value of 1.
    if (!is.null(network$learning$args$iss.mu))
      iss.mu = network$learning$args$iss.mu
    else
      iss.mu = 1

  }#ELSE

  # coerce iss to integer.
  return(as.numeric(iss.mu))

}#CHECK.ISS

# check the imaginary sample size for the Wishart prior over the covariance.
check.iss.w = function(iss.w, network, data) {

  if (!is.null(iss.w)) {

    if (!is.positive.integer(iss.w))
      stop("the imaginary sample size for the Wishart prior component must be a positive integer number.")

  }#THEN
  else {

    if (!is.null(network$learning$args$iss.w))
      iss.w = network$learning$args$iss.w
    else
      iss.w = ncol(data) + 2

  }#ELSE

  return(iss.w)

}#CHECK.ISS.W

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
      # also, convert everything to integers to make things simpler at the
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

    # warn if using a non-standard penalty.
    if (grepl("^aic", score) && (k != 1))
      warning("using AIC with a non-standard penalty k = ", k, ".")
    if (grepl("^bic", score) && (k != log(nrow(data))/2))
      warning("using BIC with a non-standard penalty k = ", k, ".")

  }#THEN
  else {

    # check whether there is a penalization coefficient stored in the bn object,
    # use the default for the score function otherwise.
    if (!is.null(network$learning$args$k) && (score == network$learning$test))
      k = network$learning$args$k
    else
      k = ifelse(grepl("^aic", score), 1, log(nrow(data))/2)

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

    # check whether prior is a string, and whether the label matches a known prior.
    check.label(prior, choices = prior.distributions, labels = prior.labels,
      argname = "prior distribution", see = score)

  }#ELSE

  return(prior)

}#CHECK.GRAPH.PRIOR

# check the sparsity parameter of the prior distribution over the graph space.
check.graph.hyperparameters = function(beta, prior, network, data,
    learning = FALSE) {

  default.beta =
    list("uniform" = NULL, "vsp" = 1/(ncol(data) - 1),
      "marginal" = structure(0.5, nodes = names(network$nodes)),
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

    if (prior %in% c("uniform")) {

      warning("unused argument beta.")
      beta = NULL

    }#THEN
    else if (prior == "marginal") {

      if (!is.probability(beta) || (beta >= 1))
       stop("beta must be a probability smaller than 1.")

      beta = structure(noattr(beta, ok = character(0)), nodes = names(network$nodes))

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

# check the number of scores to average.
check.l = function(l) {

  if (is.null(l))
    l = 5
  else
    if (!is.positive.integer(l))
      stop("l must be a positive integer, the number of scores to average.")

  return(as.numeric(l))

}#CHECK.L

# check the test data for predictive scores.
check.newdata = function(newdata, network = network, data = data) {

  if (is.null(newdata))
    stop("predictive scores require a test set passed as 'newdata'.")

  # check whether data and newdata have the same columns.
  names.data = names(data)
  names.newdata = names(newdata)

  if (length(names.data) < length(names.newdata)) {

    if (all(names.newdata %in% names(data))) {

       warning("extra variables in newdata will be ignored.")
       newdata = .data.frame.column(newdata, names.data, drop = FALSE)

    }#THEN
    else {

      stop("some variables in the data are missing in newdata.")

    }#ELSE

  }#THEN
  else if (length(names.data) > length(names.newdata)) {

    stop("some variables in the data are missing in newdata.")

  }#THEN

  # reorder the columns of newdata to match data.
  newdata = .data.frame.column(newdata, names.data, drop = FALSE)

  # check whether data and newdata have the same data types.
  types.data = lapply(data, class)
  types.newdata = lapply(newdata, class)

  for (i in seq_along(types.data))
    if (!isTRUE(all.equal(types.data[[i]], types.newdata[[i]])))
      stop("variable ", names.data[i], " has a different class in newdata.")

  # for discrete variables, also check that the levels match.
  levels.data = lapply(data, levels)
  levels.newdata = lapply(newdata, levels)

  for (i in seq_along(levels.data))
    if (any(types.data[i] %in% "factor"))
      if (!isTRUE(all.equal(levels.data[[i]], levels.newdata[[i]])))
        stop("variable ", names.data[i], " has a different levels in newdata.")

  # make sure to return a data frame with column names.
  newdata = .data.frame(newdata)
  names(newdata) = names.data

  return(newdata)

}#CHECK.NEWDATA

# check the R function implementing a custom score.
check.custom.score.function = function(fun) {

  # there is no possible default value.
  if (is.null(fun))
    stop("missing the custom score function.")

  # check the argument list.
  fun.arguments = names(formals(fun))
  if (!setequal(fun.arguments, c("node", "parents", "data", "args")) ||
      anyDuplicated(fun.arguments))
    stop("the custom score function must have signature function(node, parents, data, args).")

  return(fun)

}#CHECK.CUSTOM.SCORE.FUNCTION

# check the additional arguments' list passed to a custom score.
check.custom.score.arguments = function(args) {

  # default to an empty argument list.
  if (is.null(args))
    args = list()
  else if (!is.list(args))
    stop("the arguments for the custom score must be passed as a list.")

  return(args)

}#CHECK.CUSTOM.SCORE.ARGUMENTS
