
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
  if (score %in% c("loglik", "loglik-g"))
    return(TRUE)
  # same with AIC and BIC.
  if (score %in% c("aic", "aic-g", "bic", "bic-g"))
    return(TRUE)
  # BDe and BGe are score equivalent if they have uniform priors (e.g. BDeu and BGeu).
  else if ((score %in% c("bde", "bge")) &&
           (extra$prior %in% c("uniform", "marginal")))
    return(TRUE)

  # a conservative default.
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
  if (score %in% c("bde", "bds", "mbde", "bge"))
    extra.args$iss = check.iss(iss = extra.args$iss,
      network = network, data = data)

  # check the graph prior distribution.
  if (score %in% c("bde", "bds", "bdj", "mbde", "bdla", "bge"))
    extra.args$prior = check.graph.prior(prior = extra.args$prior,
      network = network)

  # check the sparsity parameter of the graph prior distribution.
  if (score %in% c("bde", "bds", "bdj", "mbde", "bdla", "bge"))
    extra.args$beta = check.graph.hyperparameters(beta = extra.args$beta,
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

  # check the number of scores to average.
  if (score == "bdla")
    extra.args$l = check.l(l = extra.args$l)

  check.unused.args(extra.args, score.extra.args[[score]])

  return(extra.args)

}#CHECK.SCORE.ARGS

# check the imaginary sample size.
check.iss = function(iss, network, data) {

  # check which type of data we are dealing with.
  type = data.type(data)

  if (!is.null(iss)) {

    # validate the imaginary sample size.
    if (!is.positive(iss))
      stop("the imaginary sample size must be a positive number.")
    # if iss = 1 the bge is NaN, if iss = 2 and phi = "heckerman" the
    # computation stops with the following error:
    # Error in solve.default(phi[A, A]) :
    #   Lapack routine dgesv: system is exactly singular
    if ((type == "factor") && (iss < 1))
      warning("very small imaginary sample size, the results may be affected by numeric problems.")
    else if ((type == "continuous") && (iss < 3))
      stop("the imaginary sample size must be a numeric value greater than or equal to 3.")
    else if (iss < 1)
      stop("the imaginary sample size must be greater than or equal to 1.")

  }#THEN
  else {

    # check whether there is an imaginary sample size stored in the bn object;
    # otherwise use a the de facto standard value of 10.
    if (!is.null(network$learning$args$iss))
      iss = network$learning$args$iss
    else {

      if (type == "factor")
        iss = 1
      else
        iss = 10

    }#ELSE

  }#ELSE

  # coerce iss to integer.
  return(as.numeric(iss))

}#CHECK.ISS

# check the phi defintion to be used in the bge score.
check.phi = function(phi, network, data) {

  if (!is.null(phi)) {

    check.label(phi, choices = c("heckerman", "bottcher"),
      argname = "phi definition")

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
check.graph.hyperparameters = function(beta, prior, network, data,
    learning = FALSE) {

  default.beta =
    list("uniform" = NULL, "vsp" = 1/ncol(data),
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

