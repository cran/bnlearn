
# internal label mapping functions.
.ctf = function(nodes) paste0(nodes, ".", recycle0 = TRUE)
.exo = function(nodes) paste0("u", nodes, recycle0 = TRUE)
.fac.ctf = function(nodes) sub("\\.$", "", nodes)
.fac.exo = function(nodes) sub("^u", "", nodes)

# convert a bn object into a structural causal model.
from.bn.to.scm = function(bn) {

  # special case twin bn objects, which require far less work.
  if (is(bn, "bn.twin"))
    return(from.bn.twin.to.scm.twin(bn))

  # create the node labels and check for duplicates.
  factual = names(bn$nodes)
  exogenous = .exo(factual)

  if (anyDuplicated(c(factual, exogenous)))
    stop("duplicated node labels in the 'scm' object.")

  # add the arcs from the exogenous to the factual nodes.
  arcs = matrix(c(bn$arcs[, "from"], exogenous, bn$arcs[, "to"], factual),
           ncol = 2, dimnames = list(NULL, c("from", "to")))

  # node metadata.
  node.metadata = function(node) {

    if (node %in% factual) {

      factual = character(0)
      exogenous = .exo(node)
      parents = bn$nodes[[node]]$parents
      children = bn$nodes[[node]]$children

    }#THEN
    else {

      factual = .fac.exo(node)
      exogenous = character(0)
      parents = character(0)
      children = character(0)

    }#ELSE

    list(counterfactual = character(0), exogenous = exogenous,
         factual = factual, parents = parents, children = children)

  }#NODE.METADATA

  nodes = sapply(c(factual, exogenous), node.metadata, simplify = FALSE)

  # node roles.
  roles = list(factual = factual, exogenous = exogenous,
            counterfactual = character(0))

  return(structure(list(roles = roles, nodes = nodes, arcs = arcs),
           class = "scm"))

}#FROM.BN.TO.SCM

# convert a twin network from a bn object to a structural causal model.
from.bn.twin.to.scm.twin = function(bn) {

  # extract the node roles.
  roles = bn$learning$roles

  # node metadata.
  node.metadata = function(node, roles) {

    if (node %in% roles$factual) {

      factual = character(0)
      exogenous = .exo(node)
      counterfactual = .ctf(node)
      parents = intersect(bn$nodes[[node]]$parents, roles$factual)
      children = intersect(bn$nodes[[node]]$children, roles$factual)

    }#THEN
    else if (node %in% roles$counterfactual) {

      factual = .fac.ctf(node)
      exogenous = .exo(factual)
      counterfactual = character(0)
      parents = intersect(bn$nodes[[node]]$parents, roles$counterfactual)
      children = intersect(bn$nodes[[node]]$children, roles$counterfactual)

    }#THEN
    else if (node %in% roles$exogenous) {

      factual = .fac.exo(node)
      exogenous = character(0)
      counterfactual = .ctf(factual)
      parents = character(0)
      children = character(0)

    }#THEN

    list(counterfactual = counterfactual, exogenous = exogenous,
         factual = factual, parents = parents, children = children)

  }#NODE.METADATA

  # note whether the twin network has been modified by a counterfactual.
  if (is(bn, "bn.ctf"))
    all.classes = c("scm", "scm.twin", "scm.ctf")
  else
    all.classes = c("scm", "scm.twin")
  # assemble the return value.
  nodes = sapply(names(bn$nodes), node.metadata, roles = roles, simplify = FALSE)
  twin = structure(list(roles = roles, nodes = nodes, arcs = bn$arcs),
           class = all.classes)

  return(twin)

}#FROM.BN.TWIN.TO.SCM.TWIN

# convert a structural causal model into a bn object.
from.scm.to.bn = function(scm) {

  # discarding exogenous nodes loses information for twin networks, convert
  # them to bayesian networks as-is.
  if (is(scm, "scm.twin"))
    keep.nodes = names(scm$nodes)
  else
    keep.nodes = scm$roles$factual

  # start with an empty network.
  bn = empty.graph.backend(keep.nodes)
  # carry over arcs between factual nodes (if the tail is a factual node, the
  # head must be a factual node as well).
  bn$arcs = scm$arcs[scm$arcs[, "from"] %in% keep.nodes, ]
  # update the node metadata.
  bn$nodes = cache.structure(keep.nodes, arcs = bn$arcs)

  # preserve class and node roles to make conversions easier.
  if (is(scm, "scm.twin")) {

    bn$learning$algo = "twin"
    bn$learning$roles = scm$roles
    if (is(scm, "scm.ctf"))
      class(bn) = c("bn", "bn.twin", "bn.ctf")
    else
      class(bn) = c("bn", "bn.twin")

  }#THEN

  return(bn)

}#FROM.SCM.TO.BN
