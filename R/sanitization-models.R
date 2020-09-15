
# check an object of class bn.
check.bn = function(bn) {

  if (missing(bn))
    stop("an object of class 'bn' is required.")
  if (!is(bn, "bn")) {

    stop(sprintf("%s must be an object of class 'bn'.",
           deparse(substitute(bn))))

  }#THEN

}#CHECK.BN

# check two bn's against each other.
match.bn = function(bn1, bn2) {

  # the two networks must have the same node set.
  nodes1 = names(bn1$nodes)
  nodes2 = names(bn2$nodes)

  equal = setequal(nodes1, nodes2) && (length(nodes1) == length(nodes2))

  if (!equal)
    stop("the two networks have different node sets.")

}#MATCH.BN

# check an object of class bn or bn.fit.
check.bn.or.fit = function(bn) {

  if (missing(bn))
    stop("an object of class 'bn' or 'bn.fit' is required.")
  if (!is(bn, "bn") && !is(bn, "bn.fit")) {

    stop(sprintf("%s must be an object of class 'bn' or 'bn.fit'.",
           deparse(substitute(bn))))

  }#THEN

}#CHECK.BN.OR.FIT

# check an object of class bn.fit.
check.fit = function(bn) {

  if (missing(bn))
    stop("an object of class 'bn.fit' is required.")
  if (!is(bn, "bn.fit")) {

    stop(sprintf("%s must be an object of class 'bn.fit'.",
           deparse(substitute(bn))))

  }#THEN

}#CHECK.FIT

# check the structure of a naive Bayes classifier.
check.bn.naive = function(bn) {

  # check whether it's a valid bn/bn.fit object.
  check.bn.or.fit(bn)
  # there must be a single root node, check.
  root = root.leaf.nodes(bn, leaf = FALSE)

  if (length(root) != 1)
    stop("a naive Bayes classifier can have only one root node, the training variable.")

  # cache the node labels.
  if (is(bn, "bn"))
    nodes = names(bn$nodes)
  else
    nodes = names(bn)
  # get the explanatory variables.
  explanatory = nodes[nodes != root]
  leafs = root.leaf.nodes(bn, leaf = TRUE)
  # all the explanatory variables must be leaf nodes, check.
  if (!identical(sort(explanatory), sort(leafs)))
    stop("all the explanatory variables must be leaf nodes.")
  # all the explanatory variables must have a single parent, the root node, check.
  nparents = sapply(explanatory, function(node) { length(parents(bn, node))  })

  if (any(nparents != 1))
    stop("all the explanatory variables must be children of the training variable.")

}#CHECK.BN.NAIVE

# check the structure of a naive Bayes classifier.
check.bn.tan = function(bn) {

  # check whether it's a valid bn/bn.fit object.
  check.bn.or.fit(bn)
  # there must be a single root node, check.
  root = root.leaf.nodes(bn, leaf = FALSE)

  if (length(root) != 1)
    stop("a naive Bayes classifier can have only one root node, the training variable.")

  # that root node must be the training variable, check.
  if (is(bn, "bn")) {

    # double-check just in case.
    check.nodes(bn$learning$args$training)

    nodes = names(bn$nodes)
    training = bn$learning$args$training

  }#THEN
  else {

    # double-check just in case.
    check.nodes(attr(bn, "training"))

    nodes = names(bn)
    training = attr(bn, "training")

  }#ELSE

  if (!identical(training, root))
    stop("the training node is not the only root node in the graph.")

  # get the explanatory variables.
  explanatory = nodes[nodes != root]
  # all the explanatory variables save one must have exactly two parents, check.
  nparents = sapply(explanatory, function(node) { length(parents(bn, node))  })

  if (!( (length(which(nparents == 2)) == length(explanatory) - 1) && (length(which(nparents == 1)) == 1) ))
    stop("the explanatory variables must form a tree.")

}#CHECK.BN.TAN

