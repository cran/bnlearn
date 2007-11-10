
# compute the loglikelihood of single node of a discrete
# bayesian network.
loglik.node = function(node, x, data, debug) {

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

    tab = table(data[, node])
    node.loglik = log(dmultinom(tab, ndata, tab / ndata))

  }#THEN
  # this node has at least one parent.
  else {

    # if there is only one parent, get it easy.
    if (length(node.parents) == 1) {

      config = data[, node.parents]

    }#THEN
    else {

      config = factor(apply(as.data.frame(data[,node.parents]), 
        1, paste, sep = "", collapse = ":"))

    }#ELSE

    tab = table(data[, node], config)
    joint = dmultinom(tab, prob = tab/ndata)
    cond = dmultinom(colSums(tab), prob = colSums(tab)/ndata)
    node.loglik = log(joint) - log(cond)

    if (debug) {

      cat("  > joint log-probability:", joint, "\n")
      cat("  > conditioning log-probability:", cond, "\n")

    }#THEN

  }#ELSE

  if (debug) {

    cat("  > node loglikelihood is", node.loglik, ".\n")

  }#THEN

  node.loglik

}#LOGLIK.NODE

# set an imaginary sample size to populate the a priori distribution.
imaginary.sample.size = function(x, data) {

  2 * prod(sapply(names(x$nodes),
        function(node) { nlevels(data[, node]) }))

}#IMAGINARY.SAMPLE.SIZE

# compute the dirichlet posterior density of a node.
dirichlet.node = function(node, x, data, imaginary.sample.size, debug) {

    node.params = as.vector(nparams.backend(x, data)[node])
    node.levels = nlevels(data[, node])
    node.parents = x$nodes[[node]]$parents
    parents.levels = node.params / node.levels

    # use an uniform a priori distribution.
    cprior = rep(imaginary.sample.size / node.params , node.params)
    # update the a priori distribution to the a posteriori one.
    alpha = cprior + as.vector(table(data[, c(node, node.parents)]))

    n = imaginary.sample.size 
    N = imaginary.sample.size + nrow(data)  

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
      cat("  > imaginary sample size:\n")
      print(imaginary.sample.size)

    }#THEN

    # if the node is not a root one, the a posteriori distribution depends
    # on the configuration of the parents of the node.
    if (length(node.parents) >0) {

      # if there is only one parent, get it easy.
      if (length(node.parents) == 1) {

        config = data[, node.parents]

      }#THEN
      else {

        config = factor(apply(as.data.frame(data[,node.parents]), 1, paste,
          sep = "", collapse = ":"))

      }#ELSE

      tab = colSums(table(data[, node], config))
      alphaj = rep(imaginary.sample.size / length(tab), length(tab))

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

        res = sum(lgamma(alpha) - lgamma(cprior)) + lgamma(n) - lgamma(N)

    }#ELSE

    if (debug) {

      cat("  > node dirichlet score is", res, ".\n")

    }#THEN

    res

}

