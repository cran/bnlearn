
# compute the score of a bayesian network.
network.score = function(network, data, score, extra.args, debug = FALSE) {

  sum(per.node.score(network = network, data = data, score = score,
    nodes = names(network$nodes), extra.args = extra.args,
    debug = debug))

}#NETWORK.SCORE

# compute single nodes' contributions to the network score.
per.node.score = function(network, data, score, nodes, extra.args,
    debug = FALSE) {

  if (score == "k2") {

   res = sapply(nodes, dirichlet.node, x = network, data = data,
           debug = debug)

  }#THEN
  else if (score == "bde") {

   res = sapply(nodes, dirichlet.node, x = network,
           imaginary.sample.size = extra.args$iss, data = data,
           debug = debug)

  }#THEN
  else if (score == "loglik") {

    res = sapply(nodes, loglik.node, x = network, data = data,
            debug = debug)

  }#THEN
  else if (score %in% c("aic", "bic")) {

    res = sapply(nodes, aic.node, x = network, data = data,
            k = extra.args$k, debug = debug)

  }#THEN
  else if (score == "bge") {

    res = sapply(nodes, bge.node, x = network, phi = extra.args$phi,
            data = data, imaginary.sample.size = extra.args$iss,
            debug = debug)

  }#THEN
  else if (score == "loglik-g") {

    res = sapply(nodes, gloglik.node, x = network, data = data,
            debug = debug)

  }#THEN
  else if (score %in% c("aic-g", "bic-g")) {

    res = sapply(nodes, aic.gauss.node, x = network, data = data,
            k = extra.args$k, debug = debug)

  }#THEN

  # set the names on the resulting array for easy reference.
  names(res) = nodes

  return(res)

}#PER.NODE.SCORE

# compute the loglikelihood of single node of a discrete bayesian network.
gloglik.node = function(node, x, data, debug = FALSE) {

  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* processing node", node, ".\n")

  }#THEN

  # get the parents of the node.
  node.parents = x$nodes[[node]]$parents

  if (length(node.parents) == 0) {

    res = sum(dnorm(data[, node], mean = mean(data[, node]),
                sd = sd(data[, node]), log = TRUE))

  }#THEN
  else {

    intercept = rep(1, nrow(data))
    qr.x = qr(cbind(intercept, data[, node.parents]))
    fitted = qr.fitted(qr.x, data[, node])
    sd = sd(qr.resid(qr.x, data[, node]))

    res = sum(dnorm(data[, node], mean = fitted, sd = sd, log = TRUE))

  }#ELSE

  if (debug) {

    cat("  > node loglikelihood is", res, ".\n")

  }#THEN

  return(res)

}#GLOGLIK.NODE

# compute the AIC of single node of a discrete bayesian network.
aic.gauss.node = function(node, x, data, k = 1, debug = FALSE) {

  gloglik.node(x = x, node = node, data = data, debug = debug) -
  k * nparams.gaussian.node(node = node, x = x)

}#AIC.GAUSS.NODE

# compute the loglikelihood of single node of a discrete bayesian network.
loglik.node = function(node, x, data, debug = FALSE) {

  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* processing node", node, ".\n")

  }#THEN

  # get the parents of the node.
  node.parents = x$nodes[[node]]$parents
  # cache the sample size.
  ndata = nrow(data)

  # this node is a root node.
  if (length(node.parents) == 0) {

    node.loglik = .Call("dlik",
       x = data[, node],
       lx = nlevels(data[, node]),
       length = nrow(data),
       PACKAGE = "bnlearn")

  }#THEN
  # this node has at least one parent.
  else {

    # if there is only one parent, get it easy.
    if (length(node.parents) == 1)
      config = data[, node.parents]
    else
      config = configurations(data[,node.parents])

    node.loglik = .Call("cdlik",
       x = data[, node],
       y = config,
       lx = nlevels(data[, node]),
       ly = nlevels(config),
       length = nrow(data),
       PACKAGE = "bnlearn")

  }#ELSE

  if (debug) {

    cat("  > node loglikelihood is", node.loglik, ".\n")

  }#THEN

  return(node.loglik)

}#LOGLIK.NODE

# compute the AIC of single node of a discrete bayesian network.
aic.node = function(node, x, data, k = 1, debug = FALSE) {

  loglik.node(x = x, node = node, data = data, debug = debug) -
  k * nparams.discrete.node(node = node, x = x, data = data, real = TRUE)

}#AIC.NODE

# compute the dirichlet posterior density of a node.
dirichlet.node = function(node, x, data, imaginary.sample.size = NULL, debug = FALSE) {

  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* processing node", node, ".\n")

  }#THEN

  node.parents = x$nodes[[node]]$parents

  if (length(node.parents) == 0) {

    dirichlet = .Call("dpost",
       x = data[, node],
       lx = nlevels(data[, node]),
       length = nrow(data),
       iss = imaginary.sample.size,
       debug = debug,
       PACKAGE = "bnlearn")

  }#THEN
  else {

    # if there is only one parent, get it easy.
    if (length(node.parents) == 1)
      config = data[, node.parents]
    else
      config = configurations(data[, node.parents])

    node.params = nparams.discrete.node(node, x, data, real = FALSE)

    dirichlet = .Call("cdpost",
       x = data[, node],
       y = config,
       lx = nlevels(data[, node]),
       ly = nlevels(config),
       length = nrow(data),
       iss = imaginary.sample.size,
       nparams = node.params,
       debug = debug,
       PACKAGE = "bnlearn")

  }#ELSE

  if (debug) {

    cat("  > node dirichlet score is", dirichlet, ".\n")

  }#THEN

  return(dirichlet)

}#DIRICHLET.NODE

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

