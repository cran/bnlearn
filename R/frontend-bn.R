
# get the number of parameters of the bayesian network.
nparams = function(x, data, effective = FALSE, debug = FALSE) {

  # check x's class.
  check.bn.or.fit(x)
  # check debug and  effective.
  check.logical(debug)
  check.logical(effective)

  if (is(x, "bn")) {

    # check the data are there.
    check.data(data, allow.missing = TRUE)
    # check the network against the data.
    check.bn.vs.data(x, data)
    # the number of parameters is unknown for partially directed graphs.
    if (is.pdag(x$arcs, names(x$nodes)))
      stop("the graph is only partially directed.")

    if (effective) {

      # fit the network to compute the number of non-zero parameters.
      x = bn.fit(x, data)
      return(nparams.fitted(x, effective = TRUE, debug = debug))

    }#THEN
    else {

      return(nparams.backend(x, data = data, debug = debug))

    }#ELSE

  }#THEN
  else {

    if (!missing(data))
      warning("unused argument data.")

    return(nparams.fitted(x, effective = effective, debug = debug))

  }#ELSE

}#NPARAMS

# get the number of tests/scores used in structure learning.
ntests = function(x) {

  # check x's class.
  check.bn(x)

  x$learning$ntests

}#NTESTS

# structural hamming distance.
shd = function(learned, true, wlbl = FALSE, debug = FALSE) {

  # check x's class.
  check.bn(learned)
  check.bn(true)
  # check debug and wlbl.
  check.logical(debug)
  check.logical(wlbl)
  # the two networks must have the same node set.
  match.bn(learned, true)

  structural.hamming.distance(learned = learned, true = true, wlbl = wlbl,
    debug = debug)

}#SHD

# Hamming distance.
hamming = function(learned, true, debug = FALSE) {

  # check learned's and true's class.
  check.bn(learned)
  check.bn(true)
  # check debug.
  check.logical(debug)
  # the two networks must have the same node set.
  match.bn(learned, true)

  hamming.distance(learned = learned, true = true, debug = debug)

}#HAMMING

# get the whitelist used by the learning algorithm.
whitelist = function(x) {

  # check x's class.
  check.bn(x)

  if (is.null(x$learning$whitelist))
    return(matrix(character(0), nrow = 0, ncol = 2,
      dimnames = list(NULL, c("from", "to"))))
  else
    return(x$learning$whitelist)

}#WHITELIST

# get the blacklist used by the learning algorithm.
blacklist = function(x) {

  # check x's class.
  check.bn(x)

  if (is.null(x$learning$blacklist))
    return(matrix(character(0), nrow = 0, ncol = 2,
      dimnames = list(NULL, c("from", "to"))))
  else
    return(x$learning$blacklist)

}#BLACKLIST

# reconstruct the equivalence class of a network.
cpdag = function(x, wlbl = FALSE, debug = FALSE) {

  # check x's class.
  check.bn.or.fit(x)
  # check debug and wlbl.
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

# contruct a consistent DAG extension of a PDAG.
cextend = function(x, strict = TRUE, debug = FALSE) {

  # check x's class.
  check.bn(x)
  # check debug.
  check.logical(debug)
  # check strict.
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
  if (is.pdag(cpdag$arcs, names(cpdag$nodes)))
    if (strict)
      stop("no consistent extension of ", deparse(substitute(x)), " is possible.")
    else
      warning("no consistent extension of ", deparse(substitute(x)), " is possible.")

  return(cpdag)

}#CEXTEND

# return the colliders (both shielded and unshielded) in a network.
colliders = function(x, arcs = FALSE, debug = FALSE) {

  # check x's class.
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

  # check x's class.
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

  # check x's class.
  check.bn.or.fit(x)
  check.logical(arcs)

  # go back to the network structure if needed.
  if (is(x, "bn.fit"))
    x = bn.net(x)

  colliders.backend(x = x, return.arcs = arcs, including.shielded = TRUE,
    including.unshielded = FALSE, debug = debug)

}#UNSHIELDED.COLLIDERS

# return the v-structures in a network.
vstructs = function(x, arcs = FALSE, debug = FALSE) {

  # check x's class.
  check.bn.or.fit(x)
  # check debug and arcs.
  check.logical(arcs)
  check.logical(debug)

  unshielded.colliders(x, arcs = arcs, debug = debug)

}#VSTRUCTS

# reconstruct the equivalence class of a network.
moral = function(x, debug = FALSE) {

  # check x's class.
  check.bn.or.fit(x)
  # check debug.
  check.logical(debug)

  # go back to the network structure if needed.
  if (is(x, "bn.fit"))
    x = bn.net(x)

  dag2ug.backend(x = x, moral = TRUE, debug = debug)

}#MORAL

# mutilated network used in likelihood weighting.
mutilated = function(x, evidence) {

  # check x's class.
  check.bn.or.fit(x)
  # check the evidence, disallowing non-ideal interventions if needed.
  evidence = check.evidence(evidence, graph = x,
                ideal.only = is(x, "bn.fit.gnet"))

  if (is(x, "bn"))
    return(mutilated.backend.bn(x, evidence))
  else
    return(mutilated.backend.fitted(x, evidence))

}#MUTILATED

# test d-separation.
dsep = function(bn, x, y, z) {

  # check x's class.
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
  if (!is.dag(bn$arcs, names(bn$nodes))) {

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

  # check the class of target and current.
  check.fit(target)
  check.fit(current)
  # warn about unused arguments, but silently ignore those set by
  # all.equal.list() that are not in the generic function.
  check.unused.args(list(...), c("use.names", "check.attributes"))

  equal.backend.fit(target = target, current = current, tolerance = tolerance)

}#ALL.EQUAL.BN.FIT

