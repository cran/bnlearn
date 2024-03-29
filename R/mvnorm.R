
# transform the local distributions in a Gaussian BN into the multivariate
# normal that is the global distribution.
gbn2mvnorm.backend = function(fitted){

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

  sigma = (chol %*% t(chol))[names, names, drop = FALSE]

  return(list(mu = mu, sigma = sigma))

}#GBN2MVNORM

# factorize a multivariate normal distribution into the local distributions that
# make up a Gaussian BN into the multivariate.
mvnorm2gbn.backend = function(dag, mu, sigma) {

  roots = root.leaf.nodes(dag, leaf = FALSE)
  leaves = root.leaf.nodes(dag, leaf = TRUE)
  nodes = names(dag$nodes)
  ordnodes = topological.ordering(dag)
  nnodes = length(nodes)
  params = structure(vector(nnodes, mode = "list"), names = nodes)

  # reorder both the mean vector and the covariance matrix to follow the node
  # ordering of the nodes.
  mu = mu[ordnodes, drop = FALSE]
  orig = sigma = sigma[ordnodes, ordnodes, drop = FALSE]

  # patch the covariance matrix: if a root node has variance zero, all
  # covariances are zero as well which means that we can replace the variance
  # with a positive number without changing any other local distribution.
  zero.variance.roots = names(which(diag(sigma)[roots] == 0))
  diag(sigma)[zero.variance.roots] = 1

  # patch the covariance matrix: if a leaf node has zero residual variance, the
  # covariance matrix is singular, but we can increase it without affecting any
  # other local distribution.
  diag(sigma)[leaves] = diag(sigma)[leaves] + 1

  # sort the covariance matrix using the topological ordering of the nodes,
  # take the Cholesky decomposition and transpose it to make it lower-triangular
  # like the "chol" object in gbn2mvnorm().
  chol = try(t(chol(sigma)), silent = TRUE)

  if (is(chol, "try-error")) {

    diag(sigma) = diag(sigma) + sqrt(.Machine$double.eps)
    warning("'covmat' is singular, augmenting the diagonal by ",
      format(sqrt(.Machine$double.eps), digits = 3), ".")
    chol = t(chol(sigma))

  }#THEN

  # the standard errors are on the diagonal, where we put them
  # (chol[node, node] = stderr).
  stderr = diag(chol)
  # floating point errors may have made the standard errors of leaf nodes
  # smaller than 1 by up to sqrt(.Machine$double.eps), fix that.
  stderr[leaves][stderr[leaves] < 1] = 1

  # and here we invert the quadratic form chol %*% t(chol).
  if (nnodes == 1)
    rho = matrix(0, nrow = 1, ncol = 1)
  else
    rho = diag(1, nrow = nnodes) - diag(diag(chol)) %*% solve(chol)
  dimnames(rho) = dimnames(chol)

  for (node in nodes) {

    pars = dag$nodes[[node]]$parents
    # the intercept is computed using the standard definition.
    coefs = rho[node, pars]
    intercept = mu[node] - sum(mu[pars] * coefs)
    coefs = structure(c(intercept, coefs), names = c("(Intercept)", pars))

    # restore the original variances in (zero-variance) root and leaf nodes.
    if (node %in% zero.variance.roots)
      sd = 0
    else if (node %in% leaves)
      sd = sqrt(stderr[node]^2 - 1)
    else
      sd = stderr[node]

    params[[node]] = list(coef = coefs, sd = sd)

  }#FOR

  # construct the Gaussian BN.
  fitted = custom.fit.backend(x = dag, dist = params, ordinal = character(0),
             debug = FALSE)

  return(fitted)

}#MVNORM2GBN.BACKEND

# compute a conditional multivariate distribution from the joint.
conditional.mvnorm = function(mu, sigma, to, from, value) {

  # use a pseudoinverse, computed in the same way as in the C implementation.
  decomp = svd(sigma[from, from, drop = FALSE])
  positive = decomp$d > max(.Machine$double.eps * decomp$d[1L], 0)
  if (!any(positive))
    ginv = matrix(0, nrow = length(from), ncol = length(from))
  else
    ginv = decomp$v[, positive, drop = FALSE] %*%
           ((1 / decomp$d[positive]) *
           t(decomp$u[, positive, drop = FALSE]))

  mu[to] + sigma[to, from, drop = FALSE] %*%
    ginv %*% as.numeric(value - mu[from])

}#CONDITIONAL.MVNORM

