
# the generic as method for class bn.
as.bn = function(x, debug = FALSE) {

  UseMethod("as.bn")

}#AS.BN

# get the number of paraters of the bayesian network.
nparams = function(x, data, debug = FALSE) {

  # check x's class.
  check.bn.or.fit(x)
  # check debug.
  check.logical(debug)

  if (class(x) == "bn") {

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

# structural hamming distance.
shd = function(learned, true, debug = FALSE) {

  # check x's class.
  check.bn(learned)
  check.bn(true)
  # check debug.
  check.logical(debug)
  # the two networks must have the same node set.
  nodes = names(learned$nodes);
  check.nodes(nodes, graph = true, min.nodes = length(nodes))

  structural.hamming.distance(learned = learned, true = true, debug = debug)

}#SHD

# reconstruct the equivalence class of a network.
cpdag = function(x, debug = FALSE) {

  # check x's class.
  check.bn(x)
  # check debug.
  check.logical(debug)

  cpdag.backend(x = x, debug = debug)

}#CPDAG

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

}#CPDAG

