
# compute the logliklihood score with an optional penalty.
loglik.score = function(x, data, penalized = 0, debug = FALSE) {

  logscore = sum(sapply(names(x$nodes),
    loglik.node, x = x, data = data, debug = debug))

  if (penalized == 0)
    logscore
  else
    logscore - penalized *
      sum(nparams.backend(x, data, real = TRUE))

}#LOGLIK.SCORE

# compute the loglikelihood of single node of a discrete
# bayesian network.
loglik.node = function(node, x, data, debug = FALSE) {

  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* processing node", node, ".\n")

  }#THEN

  # get the parents of the node.
  node.parents = x$nodes[[node]]$parents
  # cache the sample size.
  ndata = nrow(data)
  # initializze the loglikelihood score.
  node.loglik = 0

  # this node is a root node.
  if (length(node.parents) == 0) {

    node.loglik = .C("dlik",
       x = data[, node],
       lx = nlevels(data[, node]),
       length = nrow(data),
       result = node.loglik,
       PACKAGE = "bnlearn")$result

  }#THEN
  # this node has at least one parent.
  else {

    # if there is only one parent, get it easy.
    if (length(node.parents) == 1)
      config = data[, node.parents]
    else
      config = configurations(data[,node.parents])

    node.loglik = .C("cdlik",
       x = data[, node],
       y = config,
       lx = nlevels(data[, node]),
       ly = nlevels(config),
       length = nrow(data),
       result = node.loglik,
       PACKAGE = "bnlearn")$result

  }#ELSE

  if (debug) {

    cat("  > node loglikelihood is", node.loglik, ".\n")

  }#THEN

  node.loglik

}#LOGLIK.NODE

# compute the AIC of single node of a discrete bayesian network.
aic.node = function(node, x, data, k = 1, debug = FALSE) {

  loglik.node(x = x, node = node, data = data, debug = debug) -
  k * nparams.node(node = node, x = x, data = data, real = TRUE)

}#AIC.NODE

# compute the BDe Dirichlet score.
bde.score = function(x, data, iss = NULL, debug = FALSE) {

  if (debug)
    cat("* imaginary sample size set to", iss, ".\n")

  sum(sapply(names(x$nodes), dirichlet.node, x = x, data = data,
    imaginary.sample.size = iss, debug = debug))

}#BDE.SCORE

# compute the K2 score.
k2.score = function(x, data, debug = FALSE) {

  sum(sapply(names(x$nodes), dirichlet.node, x = x, data = data, debug = debug))

}#K2.SCORE

# compute the dirichlet posterior density of a node.
# Reverse-engineered from the deal package, with copyright:
# Copyright (C) 2002  Susanne Gammelgaard Bøttcher, Claus Dethlefsen
# Original code licenced under GPLv2 or later version.
dirichlet.node = function(node, x, data, imaginary.sample.size = NULL, debug = FALSE) {

  node.params = nparams.node(node, x, data, real = FALSE)
  node.levels = nlevels(data[, node])
  node.parents = x$nodes[[node]]$parents
  parents.levels = node.params / node.levels

  # use an uniform a priori distribution.
  if (!is.null(imaginary.sample.size)) {

    # this is for the bde() score.
    cprior = rep(imaginary.sample.size / node.params, node.params)

    n = imaginary.sample.size
    N = imaginary.sample.size + nrow(data)

  }#THEN
  else {

    # this is for the k2() score.
    cprior = rep(1, node.params)

    n = node.params
    N = node.params + nrow(data)

  }#THEN

  # update the a priori distribution to the a posteriori one.
  alpha = cprior + as.vector(table(data[, c(node, node.parents)]))

  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* processing node", node, ".\n")
    cat("  > prior distribution is:", cprior[1], "( x", length(cprior),")\n")
    cat("  > sample distribution is:\n")
    print(as.vector(table(data[, node])))
    cat("  > posterior distribution is: ( x",  length(alpha),")\n")
    print(alpha)
    cat("  > real sample size:\n")
    print(nrow(data))

  }#THEN

  # if the node is not a root one, the a posteriori distribution depends
  # on the configuration of the parents of the node.
  if (length(node.parents) > 0) {

    # if there is only one parent, get it easy.
    if (length(node.parents) == 1)
      config = data[, node.parents]
    else
      config = configurations(data[,node.parents])

    tab = colSums(table(data[, node], config))
    alphaj = rep(n / length(tab), length(tab))

    if (debug) {

      cat("  > conditional prior distribution is:", alphaj[1],
          "( x", length(alphaj) ,")\n")
      cat("  > conditional sample distribution is: ( x ", length(tab), ")\n")
      print(tab)
      cat("  > conditional posterior distribution is:\n")
      print(tab + alphaj)

    }#THEN

    res = - sum(lgamma(alphaj + tab) - lgamma(alphaj)) +
              sum(lgamma(alpha) - lgamma(cprior))

  }#THEN
  else {

      # if the node has no parents, use a fake one with a single configuration
      # to make sense of the formula of the a posteriori distribution.
      res = sum(lgamma(alpha) - lgamma(cprior)) + lgamma(n) - lgamma(N)

  }#ELSE

  if (debug) {

    cat("  > node dirichlet score is", res, ".\n")

  }#THEN

res

}#DIRICHLET.NODE

# compute the BGe Gaussian score.
bge.score = function(x, data, iss = NULL, phi = "heckerman", debug = FALSE) {

  if (debug) {

    cat("* imaginary sample size set to", iss, ".\n")
    cat("* using phi definition from", phi, ".\n")

  }#THEN

  sum(sapply(names(x$nodes), bge.node, x = x, data = data,
    imaginary.sample.size = iss, phi = phi, debug = debug))

}#BGE.SCORE

# compute the posterior density of a gaussian node.
# Reverse-engineered from the deal package, with copyright:
# Copyright (C) 2002  Susanne Gammelgaard Bøttcher, Claus Dethlefsen
# Original code licenced under GPLv2 or later version.
bge.node = function(node, x, data, imaginary.sample.size, phi = "heckerman", 
    debug = FALSE) {

  # cache the sample size.
  n = nrow(data)
  # cache the phi multiplier.
  if (phi == "bottcher")
    phi.coef = (n - 1) / n * (imaginary.sample.size - 1)
  else if (phi == "heckerman")
    phi.coef = (n - 1) / n * (imaginary.sample.size) / (imaginary.sample.size + 1) * (imaginary.sample.size - 2)
  # initialize the result.
  result = 0

  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* processing node", node, ".\n")

  }#THEN

  if (length(x$nodes[[node]]$parents) == 0) {

    if (debug) {

      cat("  > node", node, "has no parents.\n")
      cat("  > mu:", as.double(mean(data[, node])), "\n")
      cat("  > tau:", as.double(imaginary.sample.size), "\n")
      cat("  > rho:", as.double(imaginary.sample.size), "\n")
      cat("  > phi:", as.double(var(data[, node]) * phi.coef), "\n")

    }#THEN

    # Posterior for continuous node with no parents
    result = .C("postc0",
       mu = as.double(mean(data[, node])),
       tau = as.double(imaginary.sample.size),
       rho = as.double(imaginary.sample.size),
       phi = as.double(var(data[, node]) * phi.coef),
       loglik = as.double(result),
       y = as.double(data[, node]),
       n = as.integer(n),
       PACKAGE = "bnlearn")$loglik

  }#THEN
  else {

    # Posterior for continuous node with continuous parents
    node.parents = x$nodes[[node]]$parents
    n.parents = length(node.parents)

    tau.build = function(A, mu, nu, rho, phi) {

      B = A 

      A = setdiff(1:ncol(phi),B)
      if (length(A) < 1) A = TRUE

      rho.BlA = rho + length(A)
      phi.AA.inv = solve(phi[A,A])
      phi.tmp = phi[B,A] %*% phi.AA.inv
      phi.BlA = phi[B,B] - phi.tmp%*%phi[A,B]
      mu.BlA  = c(mu[B] - phi.tmp%*%mu[A], phi.tmp)
      tau.BlA.inv.11 = 1/nu + t(mu[A])%*%phi.AA.inv%*%mu[A]
      tau.BlA.inv.22 = phi.AA.inv
      tau.BlA.inv.12 = -t(mu[A] %*% phi.AA.inv)

      tau.inv = rbind(cbind(tau.BlA.inv.11, t(tau.BlA.inv.12)),
                       cbind(tau.BlA.inv.12, tau.BlA.inv.22))

      solve(tau.inv)

    }#TAU.BUILD

    tau = tau.build(A = which(node.parents %in% names(x$nodes)), 
            mu = mean(data[, node.parents]),
            nu = imaginary.sample.size,
            rho = imaginary.sample.size,
            phi = var(data[, node.parents, drop = FALSE]) * phi.coef)

    if (debug) {

      cat("  > node", node, "parents '", node.parents,"'.\n")
      cat("  > mu:", as.double(c(mean(data[, node]), rep(0, n.parents))), "\n")
      cat("  > tau:\n")
      print(tau)
      cat("  > rho:", as.double(imaginary.sample.size + n.parents), "\n")
      cat("  > phi:", as.double(var(data[, node]) * phi.coef), "\n")

    }#THEN
    result = .C("postc",
        mu = as.double(c(mean(data[, node]), rep(0, n.parents))),
        tau = as.double(t(tau)),
        rho = as.double(imaginary.sample.size + n.parents),
        phi = as.double(var(data[, node]) * phi.coef),
        loglik = as.double(result),
        y = as.double(data[, node]),
        z = as.double(t(cbind("1" = rep(1, n), data[, node.parents]))),
        n = as.integer(n),
        d = as.integer(n.parents + 1),
        PACKAGE = "bnlearn")$loglik
 
  }#THEN

  if (debug) {

    cat("  > node gaussian score is", result, ".\n")

  }#THEN

result

}#BGE.NODE

