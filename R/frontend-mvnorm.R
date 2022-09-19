
# transform the local distributions in a Gaussian BN into the multivariate
# normal that is the global distribution.
gbn2mvnorm = function(fitted) {

  # the network must be a Gaussian BN.
  check.fit(fitted)
  if (!is(fitted, "bn.fit.gnet"))
    stop("'fitted' must be an object of class 'bn.fit.gnet'.")

  gbn2mvnorm.backend(fitted)

}#GBN2MVNORM

# factorize a multivariate normal distribution into the local distributions that
# make up a Gaussian BN into the multivariate.
mvnorm2gbn = function(dag, mu, sigma) {

  # the network must be completely directed and acyclic.
  check.bn(dag)
  if (is.pdag(dag$arcs, names(dag$nodes)))
    stop("'dag' is only partially directed.")
  if (!is.acyclic(dag$arcs, names(dag$nodes), directed = TRUE))
    stop("'dag' contains cycles.")

  nodes = names(dag$nodes)
  ordnodes = topological.ordering(dag)
  nnodes = length(nodes)

  # check the covariance matrix.
  if (missing(sigma))
    stop("'sigma' is missing.")
  if (!is(sigma, "matrix") || (ncol(sigma) != nrow(sigma)) ||
       (length(dim(sigma)) != 2))
    stop("'sigma' must be a 2-dimensional square matrix.")
  if ((nrow(sigma) != nnodes) || (ncol(sigma) != nnodes))
    stop("'sigma' must be a square matrix of size ", nnodes, ".")
  if (!is.real.vector(sigma))
    stop("the elements of 'sigma' must be real numbers.")
  check.covariance(sigma)

  # the nodes in the dag must match the variables in the covariance matrix.
  dims = dimnames(sigma)
  if (is.null(dims)) {

    dimnames(sigma) = list(nodes, nodes)

  }#THEN
  else {

    if (!setequal(dims[[1]], dims[[2]]))
      stop("the row names and the column names of 'sigma' are different.")
    if (!setequal(dims[[1]], nodes))
      stop("the row names of 'sigma' do not match with the nodes in 'dag'.")
    if (!setequal(dims[[2]], nodes))
      stop("the column names of 'sigma' do not match with the nodes in 'dag'.")

  }#ELSE

  # check the mean vector
  if (missing(mu)) {

    mu = structure(rep(0, nnodes), names = ordnodes)
    warning("'mu' is missing, all variables are assumed to the centered.")

  }#THEN
  else {

    if (!is.real.vector(mu))
      stop("'mu' must be a vector of real numbers.")
    if (length(mu) != nnodes)
      stop("'mu' has ", length(mu), " elements, but 'sigma' has ", nnodes, ".")

    if (is.null(names(mu)))
      names(mu) = nodes
    else if (!setequal(names(mu), nodes))
      stop("the names of the elements of 'mu' do not match with the nodes in 'dag'.")

  }#THEN

  mvnorm2gbn.backend(dag = dag, mu = mu, sigma = sigma)

}#MVNORM2GBN

