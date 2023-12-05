
project.distributions = function(from, onto) {

  # check that the two networks span the same variables.
  if (!setequal(names(from), names(onto)))
    stop("the two bn.fit objects have different node sets.")

  if (is(onto, c("bn.fit.dnet", "bn.fit.onet", "bn.fit.donet"))) {

    # check whether gRain is loaded.
    check.and.load.package("gRain")

    # create the junction tree from the first network, then perform one query
    # per node to get the conditional distribution of each node given the
    # parents it has in the second networks.
    jtree = from.bn.fit.to.grain(from)
    onto.dag = bn.net(onto)
    local.dists =
      structure(vector(length(from), mode = "list"), names = names(from))
    for (node in names(from)) {

      parents = onto.dag$nodes[[node]]$parents
      local.dists[[node]] =
        grain.query(jtree, nodes = c(node, parents), type = "conditional")

    }#FOR

    projection = custom.fit(onto.dag, local.dists)

  }#THEN
  else if (is(onto, "bn.fit.gnet")) {

    # for Gaussian networks, get the global distribution of the first network
    # and re-factorize following the structure of the second network.
    joint = gbn2mvnorm.backend(from)
    projection = mvnorm2gbn.backend(bn.net(onto), joint$mu, joint$sigma)

  }#THEN

  return(projection)

}#PROJECT.DISTRIBUTIONS

# Kullback-Leibler divergence between a network Q and a reference network P.
kullback.leibler = function(P, Q) {

  if (is(P, c("bn.fit.dnet", "bn.fit.onet", "bn.fit.donet"))) {

    kullback.leibler.discrete(P, Q)

  }#THEN
  else if (is(P, "bn.fit.gnet")) {

    kullback.leibler.gaussian(P, Q)

  }#THEN
  else if (is(P, "bn.fit.cgnet")) {

    stop("conditional Gaussian networks are not supported yet.")

  }#THEN

}#KULLBACK.LEIBLER

# Kullback-Leibler divergence for discrete networks.
kullback.leibler.discrete = function(P, Q) {

  jtree = from.bn.fit.to.grain(P)

  # the divergence decomposes cleanly into one component for each local
  # distribution.
  sum(sapply(nodes(P), function(node) {

    pars = P[[node]]$parents
    probP = P[[node]]$prob
    probQ = Q[[node]]$prob

    if (anyNA(probP) || anyNA(probQ))
      return(NA)

    if (length(pars) > 0) {

      parprobP = grain.query(jtree, nodes = pars, type = "joint")
      dim(probP) = dim(probQ) = c(nrow(probP), length(probP) / nrow(probP))
      kl.node = 0

      for (i in seq(ncol(probP)))
        kl.node = kl.node +
          parprobP[i] * sum(probP[, i] * log(probP[, i] / probQ[, i]), na.rm = TRUE)

    }#THEN
    else {

      kl.node = sum(probP * log(probP / probQ), na.rm = TRUE)

    }#ELSE

    return(kl.node)

  }))

}#KULLBACK.LEIBLER.DISCRETE

# Kullback-Leibler divergence for Gaussian networks.
kullback.leibler.gaussian = function(P, Q) {

  # the divdergence can only be computed from the global distribution in the
  # general case: using local distributions would require matching sets of
  # fitted values in for each node in the two networks, which must have been
  # fitted from exactly the same data set (not even constructed by experts).
  mvnP = gbn2mvnorm.backend(P)
  mvnQ = gbn2mvnorm.backend(Q)

  # even if the distributions are identical, floating point errors in matrix
  # inversions will make it so that the result is not exactly zero: it is
  # cheaper to check than to zap small numbers half-way through the computation.
  if (isTRUE(all.equal(mvnP, mvnQ)))
    return(0)

  # a matrix with missing values cannot be meaningfully (pseudo)inverted, and
  # its determinant is also NA, just propagate the NA.
  if (anyNA(mvnQ$sigma) || anyNA(mvnQ$sigma))
    return(NA)

  svdQ = svd(mvnQ$sigma)
  det.sigmaQ = prod(svdQ$d)
  pos = svdQ$d > svdQ$d[1] * .Machine$double.eps
  inv.sigmaQ = svdQ$v[, pos, drop = FALSE] %*% diag(1 / svdQ$d[pos]) %*%
                 t(svdQ$u[, pos, drop = FALSE])

  0.5 * as.numeric(log(det(mvnQ$sigma) / det(mvnP$sigma)) +
    sum(diag(inv.sigmaQ %*% mvnP$sigma)) - nrow(mvnQ$sigma) +
    t(mvnP$mu - mvnQ$mu) %*% inv.sigmaQ %*% (mvnP$mu - mvnQ$mu)
  )

}#KULLBACK.LEIBLER.GAUSSIAN

