
# better implementation of has.path().
has.path = function(from, to, nodes, amat, exclude.direct = FALSE) {

  lambda0 = matrix(as.integer(nodes == from), nrow = 1)

  if (exclude.direct) amat[from, to] = amat[to, from] = 0
  if (amat[from, to] == 1) return(TRUE)

  # this cast is needed to keep lambda1 a 0/1 vector.
  lambda1 = (lambda0 %*% amat > 0) + 0

  for (jumps in seq(1, length(nodes))) {

    # this cast is needed to keep lambda2 a 0/1 vector.
    lambda2 = (lambda1 %*% amat > 0) + 0

    # if all lambda2 == 0, no more jumps are possible.
    if (all(lambda2 == 0)) return(FALSE)

    # do not allow looping because of a backward jump over 
    # an undirected arc: use lambda0 to clean up lambda2.
    lambda2 = ((lambda2 - lambda0) > 0) + 0

    if (lambda2[1, to] == 1) return(TRUE)

    lambda0 = lambda1
    lambda1 = lambda2

  }#FOR

  # return FALSE, but you should not be here
  warning("apparently the graph has a path with more steps than there are nodes.")
  return(FALSE)

}#HAS.PATH

# convert a set of neighbourhoods in an arc list.
mb2arcs = function(mb, nodes) {

  empty.mb = sapply(mb, function(x) {(length(x$nbr) == 0) || is.null(x$nbr) || identical(x$nbr, "")})
  result = do.call(rbind, lapply(nodes[!empty.mb],
               function(x) { cbind(from = x, to = mb[[x]][['nbr']]) }))

  # return an empty matrix all markov blankets are empty.
  if (is.null(result))
    matrix(1:2, ncol = 2, dimnames = list(c(""), c("from", "to")))[0,]
  else
    result

}#MB2ARCS

# get the root nodes of a network.
rootnodes.backend = function(arcs, nodes) {

  nodes[!(nodes %in% unique(arcs[, "to"]))]

}#ROOTNODES.BACKEND

# get the leaf nodes of a network.
leafnodes.backend = function(arcs, nodes) {

  nodes[!(nodes %in% unique(arcs[, "from"]))]

}#LEAFNODES.BACKEND

# get the parents of a node.
parents.backend = function(arcs, node, undirected = FALSE) {

  if (!undirected)
    arcs[(arcs[, "to"] == node) & !is.undirected(arcs), "from"]
  else
    arcs[(arcs[, "to"] == node), "from"]

}#PARENTS.BACKEND

# get the children of a node.
children.backend = function(arcs, node, undirected = FALSE) {

  if (!undirected)
    arcs[(arcs[, "from"] == node) & !is.undirected(arcs), "to"]
  else
    arcs[(arcs[, "from"] == node), "to"]

}#CHILDREN.BACKEND

# get the markov blanket of a node.
mb.backend = function(arcs, node) {

  mb = c(nbr.backend(arcs, node),
      unlist(sapply(children.backend(arcs, node),
        function(child) {

          parents.backend(arcs, node)

        }), use.names = FALSE))

  unique(mb[mb != node])

}#MB.BACKEND

# backend of nparams, the "get the number of parameters of a 
# discrete bayesian network" function. If real = TRUE this
# function returns the number of _independent_ parameters
# (on parameter of each set is set by the constraint by
# the condition \sum \pi_{i} = 1).
nparams.backend = function(x, data, real = FALSE) {

  sapply(nodes(x),
    function(node) {

      (nlevels(data[, node]) - (1 * real)) *
        prod(unlist(sapply(x$nodes[[node]]$parents,
          function(p) {

            nlevels(data[, p])

          })
        ))

    })

}#NPARAMS.BACKEND

# backend for neighbourhood detection.
nbr.backend = function(arcs, node) {

  # this includes neighbours with undirected arcs.
  unique(c(arcs[arcs[, "from"] == node, "to"], arcs[arcs[, "to"] == node, "from"]))

}#NBR.BACKEND


