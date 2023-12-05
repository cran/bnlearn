
# helper function to ensure that query results have dimensions in the same order
# as the nodes in query.
grain.query = function(jtree, nodes, type = "marginal") {

   probdist = gRain::querygrain(jtree, nodes = nodes, type = type)

   if (length(nodes) > 1) {

     if (type == "marginal") {

       # probdist is a named list.
       if (any(names(probdist) != nodes))
         probdist = probdist[nodes]

     }#THEN
     else if (type %in% c("joint", "conditional")) {

       # probdist is a multidimensional table.
       if (any(names(dimnames(probdist)) != nodes))
         probdist = aperm(probdist, perm = nodes)

     }#ELSE

   }#THEN

   return(probdist)

}#GRAIN.QUERY

# helper function to extract the childern of a node in a "grain" object.
grain.get.children = function(fitted, node) {

  names(which(sapply(fitted, function(x) { node %in% x$parents } )))

}#GRAIN.GET.CHILDREN

# convert a "bn.fit" object from bnlearn into a "grain" object.
from.bn.fit.to.grain = function(x, compile = TRUE) {

  cpt = vector(length(x), mode = "list")
  names(cpt) = names(x)

  for (node in names(x)) {

    parents = x[[node]][["parents"]]

    if (length(parents) == 0)
      f = paste("~", node)
    else
      f = paste("~", paste(c(node, parents), collapse = "+"))

    values = x[[node]][["prob"]]
    levels = dimnames(values)[[1]]

    # gRain requires completely specified CPTs, while bnlearn does not.
    if (any(is.na(values))) {

      warning("NaN conditional probabilities in ", node,
        ", replaced with a uniform distribution.")

      values[is.na(values)] = 1/dim(values)[1]

    }#THEN

    cpt[[node]] = gRain::cptable(formula(f), values = values, levels = levels)

  }#FOR

  # suppress deprecation warnings from gRbase.
  cpt.list = suppressWarnings(gRain::compileCPT(cpt))

  return(gRain::grain(cpt.list, compile = compile))

}#FROM.BN.FIT.TO.GRAIN

# convert a "grain" object from gRain into a "bn.fit" object.
from.grain.to.bn.fit = function(x) {

  # stub the bn.fit object that is the return value.
  nodes = names(x$cptlist)
  fitted = vector(length(nodes), mode = "list")
  names(fitted) = nodes

  # first pass: get parents and CPTs.
  for (node in nodes) {

    prob = structure(as.table(x[["cptlist"]][[node]]), class = "table")
    parents = names(dimnames(prob))[-1]

    # marginal tables have no dimension names in bnlearn.
    prob = cptattr(prob)

    fitted[[node]] = structure(list(node = node, parents = parents,
                       children = NULL, prob = prob), class = "bn.fit.dnode")

  }#FOR

  # third pass: get the children.
  for (node in nodes)
    fitted[[node]][["children"]] = grain.get.children(fitted, node)

  return(structure(fitted, class = c("bn.fit", determine.fitted.class(fitted))))

}#FROM.GRAIN.TO.BN.FIT

# convert a "grain" object into a "bn.fit" object after propagaing the evidence.
from.grain.to.bn.fit.with.evidence = function(x) {

  # gather the evidence that was incorporated into the network.
  ev = gRain::getEvidence(x)
  # if there is no evidence, just call the basic conversion function to preserve
  # the floating-point numeric representation of the conditional probabilities.
  if (is.null(ev))
    return(from.grain.to.bn.fit(x))

  # absorb the evidence into the potentials, otherwise any query involving
  # evidence nodes just fails silently.
  x = gRain::absorbEvidence(x)

  # stub the bn.fit object that is the return value.
  nodes = names(x$cptlist)
  fitted = vector(length(nodes), mode = "list")
  names(fitted) = nodes

  # build the conditional probability tables using gRain.
  for (node in nodes) {

    prob = structure(as.table(x[["cptlist"]][[node]]), class = "table")
    parents = names(dimnames(prob))[-1]

    cpt.dimensions = c(node, parents)

    if (any(cpt.dimensions %in% ev$nodes)) {

      if (node %in% ev$nodes) {

        prob[] = ev$evi_weight[[which(node == ev$nodes)]]
        propagated = prob

      }#THEN
      else {

        propagated =
          grain.query(x, nodes = cpt.dimensions, type = "conditional")

      }#ELSE

    }#THEN
    else {

      if (length(parents) == 0) {

        propagated =
          cptattr(grain.query(x, nodes = node, type = "marginal")[[node]])

      }#THEN
      else {

        propagated =
          grain.query(x, nodes = cpt.dimensions, type = "conditional")

      }#ELSE

    }#ELSE

    fitted[[node]] = structure(list(node = node, parents = parents,
                       children = NULL, prob = propagated),
                       class = "bn.fit.dnode")

  }#FOR

  # third pass: get the children.
  for (node in nodes)
    fitted[[node]][["children"]] = grain.get.children(fitted, node)

  return(structure(fitted, class = c("bn.fit", determine.fitted.class(fitted))))

}#FROM.GRAIN.TO.BN.FIT.WITH.EVIDENCE

# determine the amount of memory used by a gRain junction tree.
tree.size = function(jtree, node.nlevels) {

  sum(sapply(jtree$rip$cliques, function(x) prod(node.nlevels[x])))

}#TREE.SIZE

