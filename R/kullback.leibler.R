
# Shannon's entropy of a fitted network P.
shannon.entropy = function(P) {

  if (is(P, c("bn.fit.dnet", "bn.fit.onet", "bn.fit.donet"))) {

    # check whether gRain is loaded.
    check.and.load.package("gRain")

    shannon.entropy.discrete(P)

  }#THEN
  else if (is(P, "bn.fit.gnet")) {

    shannon.entropy.gaussian(P)

  }#THEN
  else if (is(P, "bn.fit.cgnet")) {

    # check whether gRain is loaded.
    check.and.load.package("gRain")

    shannon.entropy.conditional.gaussian(P)

  }#THEN

}#SHANNON.ENTROPY

# Shannon's entropy for discrete networks.
shannon.entropy.discrete = function(P) {

  entropy = structure(numeric(length(P)), names = names(P))

  # we only need exact inference for non-root nodes, so there is not point in
  # constructing the junction tree for empty networks.
  if (narcs.backend(P) > 0)
    jtree = from.bn.fit.to.grain(P)

  for (node in names(P)) {

    parents = P[[node]]$parents
    cpt = P[[node]]$prob

    if (anyNA(cpt)) {

      # the node is not completely representable, we cannot compute its entropy.
      entropy[node] = NA

    }#THEN
    else if (length(parents) == 0) {

      # the entropy of a root node follows the textbook formula.
      entropy[node] = - sum(cpt[cpt != 0] * log(cpt[cpt != 0]))

    }#THEN
    else {

      # compute the probabilities of all parents configurations...
      configuration.probs = grain.query(jtree, nodes = parents, type = "joint")
      # ... compute the entropy for each column of the CPT...
      p.log.p = cpt * log(cpt)
      p.log.p[cpt == 0] = 0
      cond.entropies = - colSums(p.log.p)
      # ...and combine them to obtain the node entropy.
      entropy[node] = sum(configuration.probs * cond.entropies)

    }#ELSE

  }#FOR

  return(sum(entropy))

}#SHANNON.ENTROPY.DISCRETE

# Shannon's entropy for Gaussian networks.
shannon.entropy.gaussian = function(P) {

  sigmas = sapply(P, `[[`, "sd")

  # the entropy of each node is a simple function of its variance.
  entropy = 0.5 * log(2 * pi * sigmas^2) + 0.5

  return(sum(entropy))

}#SHANNON.ENTROPY.GAUSSIAN

# Shannon's entropy for conditional Gaussian networks.
shannon.entropy.conditional.gaussian = function(P) {

  entropy = structure(numeric(length(P)), names = names(P))

  # separate the discrete part of the network and produce its junction tree.
  node.types = sapply(P, class)
  discrete.P = subgraph(P, names(which(node.types == "bn.fit.dnode")))

  extract.and.fix = function(x) {

    cpt = P[[x]]$prob
    # paper over unidentifiable parameters with uniform distributions, which is
    # what gRain would do in any case later on.
    cpt[is.na(cpt)] = 1 / nrow(cpt)
    return(cpt)

  }#EXTRACT.AND.FIX

  dists = sapply(nodes(discrete.P), extract.and.fix)
  discrete.P = custom.fit.backend(discrete.P, dist = dists, ordinal = FALSE)
  jtree = from.bn.fit.to.grain(discrete.P)

  for (node in names(P)) {

    if (node.types[node] == "bn.fit.dnode") {

      parents = P[[node]]$parents
      cpt = P[[node]]$prob

      if (anyNA(cpt)) {

        # the node is not completely representable, we cannot compute its entropy.
        entropy[node] = NA

      }#THEN
      if (length(parents) == 0) {

        # the entropy of a root node follows the textbook formula.
        entropy[node] = - sum(cpt[cpt != 0] * log(cpt[cpt != 0]))

      }#THEN
      else {

        # compute the probabilities of all parents configurations...
        configuration.probs = grain.query(jtree, nodes = parents, type = "joint")
        # ... compute the entropy for each column of the CPT...
        p.log.p = cpt * log(cpt)
        p.log.p[cpt == 0] = 0
        cond.entropies = - colSums(p.log.p)
        # ...and combine them to obtain the node entropy.
        entropy[node] = sum(configuration.probs * cond.entropies)

      }#ELSE

    }#THEN
    else if (node.types[node] == "bn.fit.gnode") {

      # this is a Gaussian node, it has no discrete parents: the entropy is just
      # a function of its variance.
      entropy[node] = 0.5 * log(2 * pi * P[[node]]$sd^2) + 0.5

    }#THEN
    else if (node.types[node] == "bn.fit.cgnode") {

      # extract the discrete parents...
      discrete.parents = P[[node]]$parents[P[[node]]$dparents]
      # ... compute the probabilities of their configurations...
      configuration.probs =
        grain.query(jtree, nodes = discrete.parents, type = "joint")
      # ... and combine them with the variances to obtain the node entropy.
      entropy[node] =
        sum(configuration.probs * 0.5 * log(2 * pi * P[[node]]$sd^2)) + 0.5

    }#THEN

  }#FOR

  return(sum(entropy))

}#SHANNON.ENTROPY.CONDITIONAL.GAUSSIAN

# Kullback-Leibler divergence between a network Q and a reference network P.
kullback.leibler = function(P, Q) {

  if (is(P, c("bn.fit.dnet", "bn.fit.onet", "bn.fit.donet"))) {

    kullback.leibler.discrete(P, Q)

  }#THEN
  else if (is(P, "bn.fit.gnet")) {

    kullback.leibler.gaussian(P, Q)

  }#THEN
  else if (is(P, "bn.fit.cgnet")) {

    kullback.leibler.conditional.gaussian(P, Q)

  }#THEN

}#KULLBACK.LEIBLER

# Kullback-Leibler divergence for discrete networks.
kullback.leibler.discrete = function(P, Q) {

  # the Kullback-Leibler divergence is computed as the combination of two terms:
  # LL(P, P), which is essentially the entropy of P, and LL(P, Q), which is
  # essentially the cross-entropy between P and Q.
  llpp = llpq = 0

  jtreeP = from.bn.fit.to.grain(P)

  # for LL(p, q): for each node...
  for (node in names(Q)) {

    # ... find out its parents in Q (if any)...
    familyQ = c(node, Q[[node]]$parents)
    # ... compute the associated terms in the summation: the logarithms of
    # conditional probabilities in the CPT of the node in Q...
    log.terms = log(Q[[node]]$prob)
    # ... and the logarithms of the corresponding joint probabilities in P (but
    # still using the parents from Q), computed by exact inference...
    cross.terms = grain.query(jtreeP, nodes = familyQ, type = "joint")
    # ... ensure that the dimensions match between them...
    if (length(familyQ) > 1)
      log.terms = aperm(log.terms, perm = familyQ)
    log.terms = do.call(`[`, c(list(log.terms), dimnames(cross.terms)))
    # ... and sum all the terms up.
    llpq = llpq + sum(cross.terms * log.terms)

  }#FOR

  # for LL(p, p): compute the entropy of P.
  # llpp = querygrain(jtreeP, nodes = names(P), type = "joint")
  # llpp = sum(llpp * log(llpp))
  llpp = -shannon.entropy.discrete(P)

  # combine LL(P, P) and LL(p, q) to form KL(P, Q).
  return(ifelse(isTRUE(all.equal(llpp, llpq)), 0, llpp - llpq))

}#KULLBACK.LEIBLER.DISCRETE

# Kullback-Leibler divergence for Gaussian networks.
kullback.leibler.gaussian = function(P, Q) {

  # taken from gbn2mvnorm().
  fitted.to.sparse = function(fitted) {

    names = names(fitted)
    ordnames = topological.ordering(fitted)
    nnodes = length(fitted)
    mu = structure(rep(0, nnodes), names = names)
    chol =
      matrix(0, nrow = nnodes, ncol = nnodes, dimnames = list(ordnames, ordnames))

    for (node in ordnames) {

      pars = fitted[[node]]$parents
      coefs = fitted[[node]]$coefficients
      stderr = fitted[[node]]$sd

      mu[node] = sum(c(1, mu[pars]) * coefs[c("(Intercept)", pars)])
      chol[node, node] = stderr
      chol[node, ] = chol[node, ] + t(coefs[pars]) %*% chol[pars, ]

    }#FOR

    return(list(mu = mu, C = chol))

  }#FITTED.TO.SPARSE

  # construct the Cholesky matrices from the network parameters...
  distP = fitted.to.sparse(P)
  distQ = fitted.to.sparse(Q)
  # ... check whether the two networks are actually different...
  if (isTRUE(all.equal(distP, distQ)))
    return(0)
  # ... and
  if (anyNA(distP$mu) || anyNA(distP$C) || anyNA(distQ$mu) || anyNA(distQ$C))
    return(NA_real_)
  # ... compute the determinants...
  detP = prod(diag(distP$C))^2
  detQ = prod(diag(distQ$C))^2
  # ... check that the two networks are not singular...
  if (detP == 0)
    return(Inf)
  if (detQ == 0)
    return(-Inf)
  # ... compute the Frobenius norm (it works only if the rows/columns of the two
  # matrices are ordered in the same way)...
  ooo = rownames(distP$C)
  tr = sum((zapsmall(solve(distQ$C)[ooo, ooo] %*% distP$C[ooo, ooo]))^2)
  # ... compute the quadratic form...
  vec = zapsmall(solve(distQ$C)[ooo, ooo]) %*% (distQ$mu[ooo] - distP$mu[ooo])
  quad = t(vec) %*% vec

  return(as.numeric(0.5 * (tr + quad - nnodes(P) + log(detQ / detP))))

}#KULLBACK.LEIBLER.GAUSSIAN

# Kullback-Leibler divergence for Gaussian networks.
kullback.leibler.conditional.gaussian = function(P, Q) {

  gbn.from.cgbn = function(bn, value) {

    continuous.nodes = names(which(sapply(bn,
                   function(x) class(x) %in% c("bn.fit.gnode", "bn.fit.cgnode"))))

    ldists = vector(length(continuous.nodes), mode = "list")
    names(ldists) = continuous.nodes

    # for each node...
    for (node in continuous.nodes) {

      if (length(bn[[node]]$dparents) == 0) {

        # ... if the node has no discrete parents, there is only one set of
        # parameters.
        ldists[[node]] = list(coef = bn[[node]]$coefficients, sd = bn[[node]]$sd)

      }#THEN
      else {

        # ... if the node has one or more discrete parents ...
        discrete.parents = bn[[node]]$parents[bn[[node]]$dparents]
        # ... there are multiple sets of parameters...
        configurations = expand.grid(bn[[node]]$dlevels)
        # ... and we need to pick the right one...
        lookup = (interaction(configurations) ==
                  interaction(value[, discrete.parents, drop = FALSE]))

        ldists[[node]] = list(coef = (bn[[node]]$coefficients)[, which(lookup)],
                              sd = (bn[[node]]$sd)[which(lookup)])

      }#THEN

    }#FOR

    subdag = subgraph(bn.net(bn), continuous.nodes)

    return(custom.fit.backend(subdag, ldists, ordinal = FALSE))

  }#GBN.FROM.CGBN

  # create the BNs spanning the discrete nodes.
  discrete.nodes = names(which(sapply(P, is, "bn.fit.dnode")))
  dagP = subgraph.backend(bn.net(P), discrete.nodes)
  distsP = sapply(discrete.nodes, function(node) P[[node]]$prob)
  dagQ = subgraph.backend(bn.net(Q), discrete.nodes)
  distsQ = sapply(discrete.nodes, function(node) Q[[node]]$prob)

  if (any(sapply(distsP, anyNA)) || any(sapply(distsQ, anyNA)))
    return(NA)

  subP = custom.fit.backend(dagP, distsP, ordinal = FALSE)
  subQ = custom.fit.backend(dagQ, distsQ, ordinal = FALSE)

  # compute the KL between the BNs spanning the discrete nodes.
  kl = kullback.leibler.discrete(subP, subQ)

  # find the discrete parents of the continuous nodes.
  DeltaP = unique(unlist(sapply(P, function(x) x$parents[x$dparents])))
  DeltaQ = unique(unlist(sapply(Q, function(x) x$parents[x$dparents])))
  Delta = union(DeltaP, DeltaQ)

  unidentifiable = sapply(setdiff(nodes(P), discrete.nodes),
                          function(x) anyNA(P[[x]]$coefficients) ||
                                      anyNA(P[[x]]$sd) ||
                                      anyNA(Q[[x]]$coefficients) ||
                                      anyNA(Q[[x]]$sd))

  if (any(unidentifiable))
    return(NA)

  if (length(Delta) == 0) {

    # continuous nodes have no discrete parents, so there is no mixture.
    subP = gbn.from.cgbn(P, NULL)
    subQ = gbn.from.cgbn(Q, NULL)

    kl = kl + kullback.leibler.gaussian(subP, subQ)

  }#THEN
  else {

    # compute the probabilities that weight the KL divergences between the
    # elements of the mixture.
    probs = grain.query(as.grain(subP), nodes = Delta, type = "joint")
    probs = cbind(expand.grid(dimnames(probs)), value = as.numeric(probs))

    for (i in seq(nrow(probs))) {

      subP = gbn.from.cgbn(P, probs[i, , drop = FALSE])
      subQ = gbn.from.cgbn(Q, probs[i, , drop = FALSE])

      kl = kl + probs[i, "value"] * kullback.leibler.gaussian(subP, subQ)

    }#FOR

  }#ELSE

  return(kl)

}#KULLBACK.LEIBLER.CONDITIONAL.GAUSSIAN
