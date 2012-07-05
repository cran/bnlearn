
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

   res = vapply(nodes, dirichlet.node, x = network, data = data,
           debug = debug, FUN.VALUE = template.numeric)

  }#THEN
  else if (score == "bde") {

   res = vapply(nodes, dirichlet.node, x = network,
           imaginary.sample.size = extra.args$iss,
           experimental = NULL, sparse = FALSE, data = data,
           debug = debug, FUN.VALUE = template.numeric)

  }#THEN
  else if (score == "mbde") {

   res = vapply(nodes, dirichlet.node, x = network,
           imaginary.sample.size = extra.args$iss,
           experimental = extra.args$exp, sparse = FALSE,
           data = data, debug = debug, FUN.VALUE = template.numeric)

  }#THEN
  else if (score == "bdes") {

   res = vapply(nodes, dirichlet.node, x = network,
           imaginary.sample.size = extra.args$iss,
           experimental = NULL, sparse = TRUE, data = data,
           debug = debug, FUN.VALUE = template.numeric)

  }#THEN
  else if (score == "loglik") {

    res = vapply(nodes, loglik.node, x = network, data = data,
            debug = debug, FUN.VALUE = template.numeric)

  }#THEN
  else if (score %in% c("aic", "bic")) {

    res = vapply(nodes, aic.node, x = network, data = data,
            k = extra.args$k, debug = debug, FUN.VALUE = template.numeric)

  }#THEN
  else if (score == "bge") {

    res = vapply(nodes, bge.node, x = network, phi = extra.args$phi,
            data = data, iss = extra.args$iss, debug = debug,
            FUN.VALUE = template.numeric)

  }#THEN
  else if (score == "loglik-g") {

    res = vapply(nodes, gloglik.node, x = network, data = data,
            debug = debug, FUN.VALUE = template.numeric)

  }#THEN
  else if (score %in% c("aic-g", "bic-g")) {

    res = vapply(nodes, aic.gauss.node, x = network, data = data,
            k = extra.args$k, debug = debug, FUN.VALUE = template.numeric)

  }#THEN

  # set the names on the resulting array for easy reference.
  names(res) = nodes

  return(res)

}#PER.NODE.SCORE

# compute the loglikelihood of single node of a continuous bayesian network.
gloglik.node = function(node, x, data, debug = FALSE) {

  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* processing node", node, ".\n")

  }#THEN

  # get the parents of the node.
  node.parents = x$nodes[[node]]$parents
  # extract the node's column from the data frame.
  datax = minimal.data.frame.column(data, node)

  if (length(node.parents) == 0) {

    res = sum(dnorm(datax, mean = mean(datax), sd = sd(datax), log = TRUE))

  }#THEN
  else {

    qr.x = qr(minimal.qr.matrix(data, node.parents))
    fitted = qr.fitted(qr.x, datax)
    sd = sd(qr.resid(qr.x, datax))

    res = sum(dnorm(datax, mean = fitted, sd = sd, log = TRUE))

  }#ELSE

  if (debug)
    cat("  > loglikelihood is", res, ".\n")

  return(res)

}#GLOGLIK.NODE

# compute the AIC of single node of a discrete bayesian network.
aic.gauss.node = function(node, x, data, k = 1, debug = FALSE) {

  lik = gloglik.node(x = x, node = node, data = data, debug = debug)
  pen = k * nparams.gaussian.node(node = node, x = x)

  if (debug) {

    cat("  > penalization is", pen, ".\n")
    cat("  > penalized loglikelihood is", lik - pen, ".\n")

  }#THEN

  return(lik - pen)

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
  # extract the node's column from the data frame.
  datax = minimal.data.frame.column(data, node)

  # this node is a root node.
  if (length(node.parents) == 0) {

    node.loglik = .Call("dlik",
       x = datax,
       lx = nlevels(datax),
       length = ndata,
       PACKAGE = "bnlearn")

  }#THEN
  # this node has at least one parent.
  else {

    # if there is only one parent, get it easy.
    if (length(node.parents) == 1)
      config = minimal.data.frame.column(data, node.parents)
    else
      config = configurations(minimal.data.frame.column(data, node.parents))

    node.loglik = .Call("cdlik",
       x = datax,
       y = config,
       lx = nlevels(datax),
       ly = nlevels(config),
       length = ndata,
       PACKAGE = "bnlearn")

  }#ELSE

  if (debug)
    cat("  > loglikelihood is", node.loglik, ".\n")

  return(node.loglik)

}#LOGLIK.NODE

# compute the AIC of single node of a discrete bayesian network.
aic.node = function(node, x, data, k = 1, debug = FALSE) {

  lik = loglik.node(x = x, node = node, data = data, debug = debug)
  pen = k * nparams.discrete.node(node = node, x = x, data = data, real = TRUE)

  if (debug) {

    cat("  > penalization is", pen, ".\n")
    cat("  > penalized loglikelihood is", lik - pen, ".\n")

  }#THEN

  return(lik - pen)

}#AIC.NODE

# compute the dirichlet posterior density of a node.
dirichlet.node = function(node, x, data, imaginary.sample.size = NULL,
    experimental = NULL, sparse = FALSE, debug = FALSE) {

  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* processing node", node, ".\n")

  }#THEN

  # get the parents of the node.
  node.parents = x$nodes[[node]]$parents
  # extract the node's column from the data frame.
  datax = minimal.data.frame.column(data, node)
  # extract the list of eperimental data.
  exp.data = minimal.data.frame.column(experimental, node)

  if (length(node.parents) == 0) {

    dirichlet = .Call("dpost",
       x = datax,
       iss = imaginary.sample.size,
       exp = exp.data,
       PACKAGE = "bnlearn")

  }#THEN
  else {

    config = configurations(minimal.data.frame.column(data, node.parents, drop = FALSE), all = !sparse)

    node.params = nparams.discrete.node(node, x, data, real = FALSE)

    dirichlet = .Call("cdpost",
       x = datax,
       y = config,
       iss = imaginary.sample.size,
       exp = exp.data,
       nparams = node.params,
       PACKAGE = "bnlearn")

  }#ELSE

  if (debug)
    cat("  > posterior density is", dirichlet, ".\n")

  return(dirichlet)

}#DIRICHLET.NODE

# compute the posterior density of a gaussian node.
bge.node = function(node, x, data, iss, phi = "heckerman", debug = FALSE) {

  # cache the sample size.
  n = nrow(data)
  # cache the phi multiplier.
  if (phi == "bottcher")
    phi.coef = (n - 1) / n * (iss - 1)
  else if (phi == "heckerman")
    phi.coef = (n - 1) / n * (iss) / (iss + 1) * (iss - 2)
  # get the parents of the node.
  node.parents = x$nodes[[node]]$parents

  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* processing node", node, ".\n")

  }#THEN

  if (length(node.parents) == 0) {

     result = .Call("wpost",
                    x = minimal.data.frame.column(data, node),
                    imaginary = iss,
                    phi.coef = phi.coef,
                    PACKAGE = "bnlearn")

  }#THEN
  else {

    result = .Call("cwpost",
                   x = minimal.data.frame.column(data, node),
                   z = minimal.data.frame.column(data, node.parents, drop = FALSE),
                   imaginary = iss,
                   phi.coef = phi.coef,
                   PACKAGE = "bnlearn")

  }#THEN

  if (debug)
    cat("  > posterior density is", result, ".\n")

  return(result)

}#BGE.NODE

