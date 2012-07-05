
# get the number of parameters of the bayesian network.
nparams = function(x, data, debug = FALSE) {

  # check x's class.
  check.bn.or.fit(x)
  # check debug.
  check.logical(debug)

  if (is(x, "bn")) {

    # check the data are there.
    check.data(data)
    # check the network against the data.
    check.bn.vs.data(x, data)
    # nparams is unknown for partially directed graphs.
    if (is.pdag(x$arcs, names(x$nodes)))
      stop("the graph is only partially directed.")

    if (is.data.discrete(data))
      return(nparams.discrete(x, data, real = TRUE, debug = debug))
    else
      return(nparams.gaussian(x, debug = debug))

  }#THEN
  else {

    return(nparams.fitted(x, debug = debug))

  }#ELSE

}#NPARAMS

# get the number of tests/scores used in structure learning.
ntests = function(x) {

  # check x's class.
  check.bn(x)

  x$learning$ntests

}#NTESTS

# structural hamming distance.
shd = function(learned, true, debug = FALSE) {

  # check x's class.
  check.bn(learned)
  check.bn(true)
  # check debug.
  check.logical(debug)
  # the two networks must have the same node set.
  match.bn(learned, true)

  structural.hamming.distance(learned = learned, true = true, debug = debug)

}#SHD

# hamming distance and related quantities.
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
cpdag = function(x, debug = FALSE) {

  # check x's class.
  check.bn(x)
  # check debug.
  check.logical(debug)
  # check whether the graph is acclic, to be sure to return a DAG.
  if (!is.acyclic(x$arcs, names(x$nodes)))
    stop("the specified network contains cycles.")

  cpdag.backend(x = x, debug = debug)

}#CPDAG

cextend = function(x, debug = FALSE) {

  # check x's class.
  check.bn(x)
  # check debug.
  check.logical(debug)
  # check whether the graph is acclic, to be sure to return a DAG.
  if (!is.acyclic(x$arcs, names(x$nodes)))
    stop("the specified network contains cycles.")

  cpdag.extension(x = x, debug = debug)

}#CEXTEND

# report v-structures in the network.
vstructs = function(x, arcs = FALSE, debug = FALSE) {

  # check x's class.
  check.bn(x)
  # check debug and arcs.
  check.logical(arcs)
  check.logical(debug)

  vstructures(x = x, arcs = arcs, debug = debug)

}#VSTRUCTS

# reconstruct the equivalence class of a network.
moral = function(x, debug = FALSE) {

  # check x's class.
  check.bn(x)
  # check debug.
  check.logical(debug)

  dag2ug.backend(x = x, moral = TRUE, debug = debug)

}#MORAL

# test d-separation.
dsep = function(bn, x, y, z) {

  # check x's class.
  check.bn.or.fit(bn)
  # check the sets of nodes.
  check.nodes(x, graph = bn, min.nodes = 1, max.nodes = 1)
  check.nodes(y, graph = bn, min.nodes = 1, max.nodes = 1)
  if (missing(z))
    z = c()
  else
    check.nodes(z, graph = bn, min.nodes = 0)

  # check whether x and y are disjoint from z.
  if (x %in% z)
    stop("x must not be one of the nodes in z.")
  if (y %in% z)
    stop("y must not be one of the nodes in z.")

  # go back to the network structure if needed.
  if (is(bn, "bn.fit"))
    bn = bn.net(bn)
  # if the graph is not directed, take it as a CPDAG and extend it.
  if (!is.dag(bn$arcs, names(bn$nodes)))
    bn = cpdag.extension(bn)

  dseparation(bn = bn, x = x, y = y, z = z)

}#DSEP
