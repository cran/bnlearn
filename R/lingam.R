
# direct LiNGAM causal discovery.
dlingam.backend = function(data, cluster = NULL, whitelist, blacklist, mi,
    maximize = "alasso", maximize.args = maximize.args, debug = FALSE) {

  # find the causal topological ordering.
  ordering = dlingam.ordering(data = data, whitelist = whitelist,
               blacklist = blacklist, mi = mi, cluster = cluster, debug = debug)

  # identify the direct causal effects starting from the node ordering.
  if (maximize == "alasso") {

    parent.sets =
      dlingam.alasso(data = data, ordering = ordering, whitelist = whitelist,
        blacklist = blacklist, gamma = maximize.args$gamma,
        lambda.min.ratio = maximize.args$lambda.min.ratio, cluster = cluster,
        pmax = maximize.args$pmax, k = maximize.args$k, debug = debug)

  }#THEN

  # create the DAG from the adjacency matrix.
  dag = empty.graph(names(data))
  arcs = alist2arcs(parent.sets, parents = TRUE)
  dag$arcs = arcs
  dag$nodes = cache.structure(names(dag$nodes), arcs = arcs)

  # save a copy of the learned node ordering.
  dag$learning$ordering = ordering
  # save the number of tests.
  dag$learning$ntests = test.counter()

  return(dag)

}#DLINGAM.BACKEND

# causal topological ordering through regression.
dlingam.ordering = function(data, mi, whitelist, blacklist, cluster = NULL,
    debug = FALSE) {

  nodes = candidates = names(data)
  ordering = character(0)
  clean = data

  # establish precedence in the node ordering from whitelists and blacklists.
  topological.precedence = function (wl, bl, nodes) {

    # create an auxiliary graph with the blacklisted arcs, reversed.
    amat.bl = arcs2amat(bl, nodes)
    amat.viable = (1L - amat.bl)
    diag(amat.viable) = 0L

    # check whether there is a path between every pair of nodes, to establish
    # which comes first in the node ordering.
    viable.paths = matrix(FALSE, nrow = length(nodes), ncol = length(nodes),
                     dimnames = list(nodes, nodes))

    for (n1 in nodes)
      for (n2 in nodes) {

        if (n1 == n2)
          next

        viable.paths[n1, n2] =
          has.path(from = n1, to = n2, nodes = nodes, amat = amat.viable,
            exclude.direct = FALSE, underlying.graph = FALSE, debug = FALSE)

      }#FOR

    # a clear precedence: a viable path in one direction but not the other.
    preceding = (viable.paths == TRUE) & (t(viable.paths) == FALSE)
    preceding = apply(preceding, 2, function(x) nodes[x], simplify = FALSE)
    # construct the list of root nodes (they have no viable incoming path).
    roots = nodes[colSums(viable.paths) == 0]
    # construct the list of leaf nodes (they have no viable outgoing path).
    leaves = nodes[rowSums(viable.paths) == 0]
    # make sure that isolated nodes, which are both roots and leaves, are only
    # listed among the root nodes.
    leaves = setdiff(leaves, roots)

    # scan the directed arcs in the whitelist.
    wl = wl[which.directed(wl, nodes), , drop = FALSE]
    for (n1 in nodes) {

      # each arc defines a parent-child relationship, and all parents must come
      # before their children in the topological ordering.
      parents = wl[wl[, "to"] == n1, "from"]
      if (length(parents) > 0)
        preceding[[n1]] = union(preceding[[n1]], parents)

    }#FOR

    return(list(preceding = preceding, roots = roots, leaves = leaves))

  }#TOPOLOGICAL.PRECEDENCE

  topo = topological.precedence(whitelist, blacklist, nodes)

  # the roots necessarily come first.
  ordering = topo$roots
  # the candidates to be ordered are the remaining non-leaf nodes.
  candidates = setdiff(candidates, c(ordering, topo$leaves))

  if (debug)
    cat("* initial node ordering is:", ordering, ".\n")

  # until the ordering is complete...
  for (f in seq(ncol(data))) {

    viable.candidates = candidates

    # a node is a viable candidate if all preceding nodes have already been
    # included in the topological ordering.
    for (node in candidates) {

      viable = all(topo$preceding[[node]] %in% ordering)
      if (viable)
        next

      viable.candidates = setdiff(viable.candidates, node)

      if (debug)
        cat("* skipping node", node, "because these are not in the ordering yet:",
          topo$preceding[[node]], ".\n")

    }#THEN

    # ... find the node to put in the next slot...
    found = search.exogenous.variable(data = clean, to.test = candidates,
              candidates = viable.candidates, mi = mi, cluster = cluster,
              debug = debug)
    # ... remove its information from the remaining candidates...
    for (i in candidates)
      if (i != found)
        clean[, i] = dlingam.remove.effect(clean[, i], clean[, found])
    # ... and update both the ordering and the candidate set.
    ordering = c(ordering, found)
    candidates = setdiff(candidates, found)

    if (length(candidates) == 0)
      break

    if (debug)
      cat("@ current causal ordering is:", ordering, ".\n")

  }#FOR

  # leaves necessarily come last.
  if (length(topo$leaves) > 0) {

    ordering = c(ordering, topo$leaves)

    if (debug)
      cat("* adding leaf nodes", topo$leaves, ".\n")

  }#THEN

  if (debug)
    cat("@ final causal ordering is:", ordering, ".\n")

  return(ordering)

}#DLINGAM.ORDERING

# identify exogenous variables by (lack of) association with other variables.
search.exogenous.variable = function(data, candidates, to.test, mi,
  cluster = NULL, debug = FALSE) {

  per.candidate = function(i, data, to.test, debug) {

     if (debug)
       cat("> considering candidate", i, "\n")

     total.mi = 0
     xi = as.numeric(.scale(data[, i]))

     # ... for each of the nodes that are not in the topological ordering yet...
     for (j in to.test) {

       if (i == j)
         next

       xj = as.numeric(.scale(data[, j]))

       if (mi == "pwling") {

         # Empirical mutual information for non-Gaussian structural equation
         # models.
         mi.contrib = min(c(0, approx.mutual.information(xi, xj)))^2

       }#THEN
       else if (mi == "gkernel") {

         # Gaussian kernel mutual information.
         mi.contrib = kernel.mutual.information(xi, xj)

       }#THEN

       increment.test.counter()

       if (debug)
         cat("  > node", j, "has mutual information", mi.contrib, ".\n")

       total.mi = total.mi + mi.contrib

     }#FOR

     return(total.mi)

  }#PER.CANDIDATE

  association = smartSapply(cluster, candidates, per.candidate, data = data,
                  to.test = to.test, debug = debug)

  return(names(which.min(association)))

}#SEARCH.EXOGENOUS.VARIABLE

# remove the effect of one variable from another.
dlingam.remove.effect = function(xi, xj) {

  # include the intercept to allow for non-centered variables.
  as.numeric(residuals(lm(xi ~ xj)))

}#DLINGAM.REMOVE.EFFECT

# kernel mutual information from Shimizu's implementation.
kernel.mutual.information = function(x1, x2) {

  n = length(x1)

  # Hardcode the kernel parameter values suggested by Shimizu.
  if (n > 1000) {

    kappa = 2e-3
    sigma = 0.5

  }#THEN
  else {

    kappa = 2e-2
    sigma = 1.0

  }#ELSE

  # np.tile()
  X1 = matrix(rep(x1, each = n), nrow = n, ncol = n)
  X2 = matrix(rep(x2, each = n), nrow = n, ncol = n)

  # this overflows too easily with 400+ observations.
  # np.exp(-1 / (2 * sigma * 2) * (X1^2 + t(X1)^2 - 2 * X1 * t(X1)))
  K1 = exp(-1 / (2 * sigma^2) * (X1^2 + t(X1)^2 - 2 * tcrossprod(X1)))
  K2 = exp(-1 / (2 * sigma^2) * (X2^2 + t(X2)^2 - 2 * tcrossprod(X2)))

  # np.identity(n)
  In = diag(n)

  tmp1 = K1 + (n * kappa * In) / 2
  tmp2 = K2 + (n * kappa * In) / 2

  # np.r_[np.c_[tmp1 @ tmp1, K1 @ K2], np.c_[K2 @ K1, tmp2 @ tmp2]]
  K1x2 = K1 %*% K2
  tmp1x1 = tmp1 %*% tmp1
  tmp2x2 = tmp2 %*% tmp2

  K.kappa = rbind(cbind(tmp1x1, K1x2), cbind(t(K1x2), tmp2x2))

  # np.r_[np.c_[tmp1 @ tmp1, np.zeros([n, n])],
  #       np.c_[np.zeros([n, n]), tmp2 @ tmp2]]
  D.kappa = rbind(cbind(tmp1x1, matrix(0, nrow = n, ncol = n)),
                  cbind(matrix(0, nrow = n, ncol = n), tmp2x2))

  # can we use eigen() instead, since it's twice as fast?
  # np.linalg.svd(..., compute_uv = False)
  sigma.K = svd(K.kappa, nu = 0, nv = 0)$d
  sigma.D = svd(D.kappa, nu = 0, nv = 0)$d

  # (-1 / 2) * (np.sum(np.log(sigma.K)) - np.sum(np.log(sigma.D)))
  return(-0.5 * (sum(log(sigma.K)) - sum(log(sigma.D))))

}#KERNEL.MUTUAL.INFORMATION

# approximate mutual information from Hyvarinen and Smith.
approx.mutual.information = function(xi, xj) {

  rij = dlingam.remove.effect(xi, xj)
  rji = dlingam.remove.effect(xj, xi)
  sdij = cgsd(rij)
  sdji = cgsd(rji)

  # if one variable is constant, mutual information is zero by construction.
  if ((sdij == 0) || (sdji == 0))
    return(0)

  entropy = function(u) {

    k1 = 79.047
    k2 = 7.4129
    gamma = 0.37457

    term1 = (1 + log(2 * pi)) / 2
    term2 = k1 * (mean(log(cosh(u))) - gamma)^2
    term3 = k2 * (mean(u * exp((-u^2) / 2)))^2

    return(term1 - term2 - term3)

  }#ENTROPY

  term1 = entropy(xj) + entropy(rij / sdij)
  term2 = entropy(xi) + entropy(rji / sdji)

  return(term1 - term2)

}#APPROX.MUTUAL.INFORMATION

# identify sparse parent sets with adaptive lasso, given topological ordering.
dlingam.alasso = function(data, ordering, whitelist, blacklist, gamma,
    lambda.min.ratio, k, pmax = ncol(data) + 1, cluster = NULL, debug = TRUE) {

  check.and.load.package("glmnet")

  sparse.parent.set =
    function(target, ordering, data, whitelist, blacklist, debug) {

    # the nodes earlier in the topological ordering are the candidate parents.
    pos = which(ordering == target)

    if (debug)
      cat("* learning the parents of node", target, "(position", pos, ").\n")

    if (pos == 1)
      return(character(0))

    candidates = ordering[1:(pos - 1)]

    # remove blacklisted nodes from the candidates.
    if (!is.null(blacklist)) {

      blacklisted.parents = blacklist[blacklist[, "to"] == target, "from"]
      candidates = setdiff(candidates, blacklisted.parents)

      if (debug)
        cat("  > blacklisted parents:", blacklisted.parents, ".\n")

      if (length(candidates) == 0)
        return(character(0))

    }#THEN

    # it is impossible to screen a single candidate with adaptive lasso, use a
    # classic linear regression and BIC instead.
    if (length(candidates) == 1) {

      base = empty.graph.backend(c(target, candidates))
      base.bic = per.node.score(base, data = data, score = "bic-g",
                   targets = target, extra.args = list(k = k), debug = FALSE)
      arc = arc.operations(base, from = candidates, to = target, op = "set",
              check.cycles = FALSE, check.illegal = FALSE, debug = FALSE)
      arc.bic = per.node.score(arc, data = data, score = "bic-g",
                  targets = target, extra.args = list(k = k), debug = FALSE)

      if (arc.bic > base.bic)
        parents = candidates
      else
        parents = character(0)

      if (debug) {

        if (length(parents) == 0)
          cat("  > linear regression found no parents.\n")
        else
          cat("  > linear regression found parent:", parents, ".\n")

      }#THEN

      return(parents)

    }#THEN

    # shortcuts.
    y = data[, target]
    x = as.matrix(data[, candidates])

    scoring = function(i, model, k) {

      if (is.na(model$lambda[i]))
        return(-Inf)

      df = as.integer(model$df[i])
      resid = residuals(model)[, i]
      logl = dnorm(resid, mean = 0, sd = cgsd(resid, p = df), log = TRUE)

      return(sum(logl) - df * k)

    }#SCORING

    # use ridge regression to estimate the regression coefficients of the
    # candidates, which are the weights of adaptive lasso.
    rr = glmnet::glmnet(y = y, x = x, alpha = 1)
    rr$residuals = as.vector(y) - glmnet::predict.glmnet(rr, newx = x, type = "response")
    scores = sapply(seq_along(rr$lambda), scoring, model = rr, k = k)
    coefs = coefficients(rr)[-1, which.max(scores)]
    weights = 1 / abs(coefs)^gamma

    # set the weights of whitelisted parents to zero to ensure inclusion.
    if (!is.null(whitelist)) {

      whitelisted.parents = whitelist[whitelist[, "to"] == target, "from"]
      weights[whitelisted.parents] = 0

      if (debug)
        cat("  > whitelisted parents:", whitelisted.parents, ".\n")

      # all candidates are whitelisted, nothing more to do.
      if (all(weights == 0))
        return(candidates)

    }#THEN

    if (all(is.infinite(weights))) {

      # infinite weights effectively blacklist nodes in glmnet.
      parents = character(0)

    }#THEN
    else {

      # estimate the adaptive lasso coefficients by maximum BIC.
      alasso = suppressWarnings(glmnet::glmnet(y = y, x = x, alpha = 1,
                 penalty.factor = weights, pmax = pmax,
                 lambda.min.ratio = lambda.min.ratio))
      alasso$residuals =
        as.vector(y) - glmnet::predict.glmnet(alasso, newx = x, type = "response")
      scores = sapply(seq_along(alasso$lambda), scoring, model = alasso, k = k)
      coefs = coefficients(alasso)[-1, which.max(scores)]
      parents = names(which(abs(coefs) > sqrt(.Machine$double.eps)))

    }#ELSE

    if (debug) {

      if (length(parents) == 0)
        cat("  > adaptive lasso found no parents.\n")
      else
        cat("  > adaptive lasso found parents:", parents, ".\n")

      cat("  > local score:", max(scores), ".\n")

    }#THEN

    # increment the test counter by the number of penalised regressions.
    increment.test.counter(2)

    return(parents)

  }#SPARSE.PARENT.SET

  smartSapply(cluster, ordering, sparse.parent.set, ordering = ordering,
    data = data, whitelist = whitelist, blacklist = blacklist, debug = debug)

}#DLINGAM.ALASSO

