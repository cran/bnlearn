
# estimate the expected value and the covariance matrix of a multivariate
# Bernoulli distribution.
mvber.moments.backend = function(data, R, m, algorithm, algorithm.args,
    arcs, debug = FALSE) {

  nodes = names(data)
  labels = apply(unique.arcs(arcs, nodes), 1, paste, collapse = "~")

  if (is.null(arcs)) {

    nnodes = length(nodes)
    joint = matrix(0, nrow = nnodes * (nnodes - 1)/2, ncol = nnodes * (nnodes - 1)/2)

  }#THEN
  else {

    narcs = nrow(arcs)
    joint = matrix(0, nrow = narcs, ncol = narcs)

  }#ELSE

  for (r in seq_len(R)) {

    if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("* bootstrap replicate", r, ".\n")

    }#THEN

    replicate = data[sample(nrow(data), m, replace = TRUE), , drop = FALSE]

    if (debug)
      cat("* learning bayesian network structure.\n")

    # learn the network structure from the bootstrap sample.
    net = do.call(algorithm, c(list(x = replicate), algorithm.args))

    # update the counters in the matrix.
    # BEWARE: in-place modification of joint.
    .Call("mvber_joint_counters",
          arcs = net$arcs,
          nodes = nodes,
          joint = joint,
          debug = debug,
          PACKAGE = "bnlearn")

  }#FOR

  # build the return value (a list of length 2 containing the expected value
  # and the covariance matrix)
  .Call("mvber_moments",
        joint = joint,
        R = R,
        labels = labels,
        PACKAGE = "bnlearn")

}#MVBER.MOMENTS.BACKEND

# compute some descriptive statistics of a network's variability.
mvber.var.backend = function(x, method) {

  res = .Call("mvber_variance",
              var = x,
              method = method,
              PACKAGE = "bnlearn")

  return(structure(res, names = c("statistic", "normalized")))

}#MVBER.VAR.BACKEND
