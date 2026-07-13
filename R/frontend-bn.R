
# dummy methods for error handling.
residuals.bn = sigma.bn = fitted.bn = coef.bn = function(object, ...) {

  stop("not applicable to 'bn' objects.")

}#RESIDUALS.BN

# get the number of parameters of the bayesian network.
nparams = function(x, data, estimator = NULL, debug = FALSE) {

  check.bn.or.fit(x)
  check.logical(debug)

  if (is(x, "bn")) {

    # check the data are there.
    data = check.data(data, allow.missing = TRUE)
    # check the network against the data.
    data = check.bn.vs.data(x, data, reorder = TRUE)
    # the number of parameters is unknown for partially directed graphs.
    if (!is.completely.directed(x))
      stop("the graph is only partially directed.")

    # sanitize the estimator.
    estimator = check.fitting.method(estimator, data)
    check.arcs.against.assumptions(x$arcs, data, estimator)

    return(nparams.backend(x, data = data, estimator = estimator,
             debug = debug))

  }#THEN
  else {

    if (!missing(data))
      warning("unused argument data.")

    return(nparams.fitted(x, debug = debug))

  }#ELSE

}#NPARAMS

# get the number of tests/scores used in structure learning.
ntests = function(x) {

  check.bn(x)

  x$learning$ntests

}#NTESTS

# compute the Structural Hamming Distance (HD).
shd = function(learned, true, wlbl = FALSE, cpdag = TRUE, debug = FALSE) {

  check.bn(learned)
  check.bn(true)
  check.logical(debug)
  check.logical(wlbl)
  check.logical(cpdag)
  # the two networks must have the same node set.
  match.bn(learned, true)

  structural.hamming.distance(learned = learned, true = true, wlbl = wlbl,
    cpdag = cpdag, debug = debug)

}#SHD

# compute the Hamming distance.
hamming = function(learned, true, debug = FALSE) {

  check.bn(learned)
  check.bn(true)
  check.logical(debug)
  # the two networks must have the same node set.
  match.bn(learned, true)

  hamming.distance(learned = learned, true = true, debug = debug)

}#HAMMING

# compute the Structural Interventional Distance (SID).
sid = function(learned, true, debug = FALSE) {

  check.bn(learned)
  check.bn(true)
  check.logical(debug)
  # the two networks must have the same node set.
  match.bn(learned, true)

  # check whether true and learned are DAGs or CPDAGs.
  true.is.dag = directed(true) && acyclic(true)
  if (true.is.dag)
    true.is.cpdag = FALSE
  else
    true.is.cpdag = valid.cpdag.backend(true)

  learned.is.dag = directed(learned) && acyclic(learned)
  if (learned.is.dag)
    learned.is.cpdag = FALSE
  else
    learned.is.cpdag = valid.cpdag.backend(learned)

  if (!true.is.dag && !true.is.cpdag)
    stop("'true' is neither a DAG nor a valid CPDAG.")
  if (!learned.is.dag && !learned.is.cpdag)
    stop("'learned' is neither a DAG nor a valid CPDAG.")

  # the nodes must be stored in the same order to ensure consistency.
  learned$nodes = learned$nodes[names(true$nodes)]

  if (true.is.dag && learned.is.dag) {

    value = sid.dag.vs.dag(learned, true, debug = debug)

  }#THEN
  else if (true.is.dag && learned.is.cpdag) {

    eq.dags = cextend.all.backend(learned)
    value = sapply(eq.dags, sid.dag.vs.dag, true = true, debug = debug)

  }#THEN
  else if (true.is.cpdag && learned.is.dag) {

    eq.dags = cextend.all.backend(true)
    value = sapply(eq.dags, sid.dag.vs.dag, learned = learned, debug = debug)

  }#THEN
  else {

    eq.true = cextend.all.backend(true)
    eq.learned = cextend.all.backend(learned)

    value = vector(length(eq.true), mode = "list")
    for (i in seq_along(eq.true))
      value[[i]] = sapply(eq.learned, sid.dag.vs.dag, true = eq.true[[i]])

    value = unlist(value)

  }#ELSE

  # store the range of SID values for compatibility.
  attr(value, "bounds") = range(value)

  return(value)

}#SID

# get the whitelist used by the learning algorithm.
whitelist = function(x) {

  check.bn(x)

  if (is.null(x$learning$whitelist))
    return(matrix(character(0), nrow = 0, ncol = 2,
      dimnames = list(NULL, c("from", "to"))))
  else
    return(x$learning$whitelist)

}#WHITELIST

# get the blacklist used by the learning algorithm.
blacklist = function(x) {

  check.bn(x)

  if (is.null(x$learning$blacklist))
    return(matrix(character(0), nrow = 0, ncol = 2,
      dimnames = list(NULL, c("from", "to"))))
  else
    return(x$learning$blacklist)

}#BLACKLIST

# reconstruct the equivalence class of a network.
cpdag = function(x, wlbl = FALSE, debug = FALSE) {

  check.bn.or.fit(x)
  check.logical(debug)
  check.logical(wlbl)

  # go back to the network structure if needed.
  if (is(x, "bn.fit"))
    x = bn.net(x)

  # check whether the graph is acyclic, to be sure to return a DAG.
  if (!is.acyclic(x$arcs, names(x$nodes), directed = TRUE))
    stop("the specified network contains cycles.")
  # check that illegal arcs are not present in the DAG.
  if (any(which.listed(x$arcs, x$learning$illegal)))
    stop("illegal arcs present in the graph.")
  # if whitelist and blacklist are to be enforced, check that they are
  # consistent with the graph.
  if (wlbl)
    if (any(which.listed(x$arcs, x$learning$blacklist)))
      stop("blacklisted arcs present in the graph.")

  cpdag.backend(x = x, wlbl = wlbl, debug = debug)

}#CPDAG

# construct a consistent DAG extension of a PDAG.
cextend = function(x, strict = TRUE, debug = FALSE) {

  check.bn(x)
  check.logical(debug)
  check.logical(strict)
  # check whether the graph is acyclic, to be sure to return a DAG.
  if (!is.acyclic(x$arcs, names(x$nodes), directed = TRUE))
    stop("the specified network contains cycles.")
  # if the graph was learned from data, check whether arc directions have been
  # learned at all; trying to extend a skeleton (instead of a CPDAG) is probably
  # not meaningful.
  if (!is.null(x$learning$undirected) && x$learning$undirected)
    warning(deparse(substitute(x)), " is just a skeleton (no arc directions ",
      "have been learned) and trying to extend it is probably wrong.")

  cpdag = cpdag.extension(x = x, debug = debug)

  # if the graph is not completely directed, the extension was not successful.
  if (!is.completely.directed(cpdag))
    if (strict)
      stop("no consistent extension of ", deparse(substitute(x)), " is possible.")
    else
      warning("no consistent extension of ", deparse(substitute(x)), " is possible.")

  return(cpdag)

}#CEXTEND

# produce all possible extensions of a CPDAG, PDAG or MPDAG.
cextend.all = function(x, debug = FALSE) {

  check.bn(x)
  check.logical(debug)

  # the graph must be acyclic to have any consistent extension.
  if (!is.acyclic(arcs = x$arcs, nodes = names(x$nodes), directed = TRUE))
    stop("the graph 'x' contains directed cycles.")

  # a valid CPDAG uses the dedicated enumerator; a PDAG or MPDAG (a graph that
  # carries background knowledge, that is, directed arcs which are not all
  # compelled) uses the bucket enumerator, which checks its own extendability.
  if (valid.cpdag.backend(x))
    cextend.all.backend(x = x, debug = debug)
  else
    cextend.all.pdag.backend(x = x, debug = debug)

}#CEXTEND.ALL

# return the colliders (both shielded and unshielded) in a network.
colliders = function(x, arcs = FALSE, debug = FALSE) {

  check.bn.or.fit(x)
  check.logical(arcs)

  # go back to the network structure if needed.
  if (is(x, "bn.fit"))
    x = bn.net(x)

  colliders.backend(x = x, return.arcs = arcs, including.shielded = TRUE,
    including.unshielded = TRUE, debug = debug)

}#COLLIDERS

# return the unshielded colliders in a network.
unshielded.colliders = function(x, arcs = FALSE, debug = FALSE) {

  check.bn.or.fit(x)
  check.logical(arcs)

  # go back to the network structure if needed.
  if (is(x, "bn.fit"))
    x = bn.net(x)

  colliders.backend(x = x, return.arcs = arcs, including.shielded = FALSE,
    including.unshielded = TRUE, debug = debug)

}#UNSHIELDED.COLLIDERS

# return the shielded colliders in a network.
shielded.colliders = function(x, arcs = FALSE, debug = FALSE) {

  check.bn.or.fit(x)
  check.logical(arcs)

  # go back to the network structure if needed.
  if (is(x, "bn.fit"))
    x = bn.net(x)

  colliders.backend(x = x, return.arcs = arcs, including.shielded = TRUE,
    including.unshielded = FALSE, debug = debug)

}#SHIELDED.COLLIDERS

# return the v-structures in a network.
vstructs = function(x, arcs = FALSE, debug = FALSE) {

  check.bn.or.fit(x)
  check.logical(arcs)
  check.logical(debug)

  unshielded.colliders(x, arcs = arcs, debug = debug)

}#VSTRUCTS

# reconstruct the equivalence class of a network.
moral = function(x, debug = FALSE) {

  check.bn.or.fit(x)
  check.logical(debug)

  # go back to the network structure if needed.
  if (is(x, "bn.fit"))
    x = bn.net(x)

  dag2ug.backend(x = x, moral = TRUE, debug = debug)

}#MORAL

# test d-separation.
dsep = function(bn, x, y, z) {

  check.bn.or.fit(bn)
  # check the sets of nodes.
  check.nodes(x, graph = bn, max.nodes = 1)
  check.nodes(y, graph = bn, max.nodes = 1)
  if (missing(z))
    z = c()
  else
    check.nodes(z, graph = bn, min.nodes = 0)

  # go back to the network structure if needed.
  if (is(bn, "bn.fit"))
    bn = bn.net(bn)
  # if the graph is not directed, take it as a CPDAG and extend it.
  if (!is.completely.directed(bn)) {

    # trying to extend a skeleton (instead of a CPDAG) is probably not
    # meaningful.
    if (!is.null(bn$learning$undirected) && bn$learning$undirected)
      warning(deparse(substitute(x)), " is just a skeleton (no arc directions ",
        "have been learned) and trying to extend it is probably wrong.")

    bn = cpdag.extension(bn)

  }#THEN

  dseparation(bn = bn, x = x, y = y, z = z)

}#DSEP

# test the equality of two fitted networks.
all.equal.bn.fit = function(target, current, ...,
    tolerance = sqrt(.Machine$double.eps)) {

  check.fit(target)
  check.fit(current)
  # warn about unused arguments, but silently ignore those set by
  # all.equal.list() that are not in the generic function.
  check.unused.args(list(...), c("use.names", "check.attributes"))

  equal.backend.fit(target = target, current = current, tolerance = tolerance)

}#ALL.EQUAL.BN.FIT

