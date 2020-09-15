
# helper function to extract the childern of a node in a "grain" object.
grain.get.children = function(fitted, node) {

  names(which(sapply(fitted, function(x) { node %in% x$parents } )))

}#GRAIN.GET.CHILDREN

# convert a "bn.fit" object from bnlearn into a "grain" object.
from.bn.fit.to.grain = function(x) {

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

  return(gRain::grain(gRain::compileCPT(cpt)))

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

  # build the conditional probability tables using querygrain().
  for (node in nodes) {

    prob = structure(as.table(x[["cptlist"]][[node]]), class = "table")
    parents = names(dimnames(prob))[-1]

    cpt.dimensions = c(node, parents)

    if (any(cpt.dimensions %in% ev$nodes)) {

      if (node %in% ev$nodes) {

        prob[] = ev$evidence[[which(node == ev$nodes)]]
        propagated = prob

      }#THEN
      else {

        propagated = gRain::querygrain(x, nodes = cpt.dimensions, type = "conditional")

      }#ELSE

    }#THEN
    else {

      if (length(parents) == 0)
        propagated = cptattr(gRain::querygrain(x, nodes = node)[[node]])
      else
        propagated = gRain::querygrain(x, nodes = cpt.dimensions, type = "conditional")

    }#ELSE

    # make sure the dimensions of the conditional probability table are in the
    # right order, as gRain seems to shuffle them at random.
    if (length(dim(propagated)) > 1)
      propagated = aperm(propagated, match(cpt.dimensions, names(dimnames(propagated))))

    fitted[[node]] = structure(list(node = node, parents = parents,
                       children = NULL, prob = propagated),
                       class = "bn.fit.dnode")

  }#FOR

  # third pass: get the children.
  for (node in nodes)
    fitted[[node]][["children"]] = grain.get.children(fitted, node)

  return(structure(fitted, class = c("bn.fit", determine.fitted.class(fitted))))

}#FROM.GRAIN.TO.BN.FIT.WITH.EVIDENCE
