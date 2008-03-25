
# a better implementation of has.path(); it implements a one-step
# backtracking, so while it maybe perfectly fine for directed graphs
# it may report false positives (_never_ half negatives) for
# partially directed graphs.
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

  # return FALSE, but you should not be here at all.
  warning("apparently the graph has a path with more steps than there are nodes.")
  return(FALSE)

}#HAS.PATH

# a has.path() implementation for partially directed graphs.
# lots of code straight from how.many.loops().
has.pdag.path = function(from, to, nodes, amat, exclude.direct = FALSE,
    debug = FALSE) {

  if (debug) {

    cat("----------------------------------------------------------------\n")
    cat("* path discovery for arc", from, "->", to, ".\n")

  }#THEN

  # honour the exclude.direct parameter.
  if (exclude.direct) amat[from, to] = amat[to, from] = 0

  # first try the dirty way.
  if (!has.path(from, to, nodes, amat)) {

    if (debug) cat("  > no path from", from, "to", to, ".\n");

    return(FALSE)

  }#THEN

  # initialize the main buffer of the paths ...
  buffer = list(from)

  # ... a token to be set from within the sapply() call ...
  path.found = FALSE

  # ... and loop (up to #nodes - 1 iterations).
  for (i in 1:(length(nodes) - 1)) {

    # set up a temporary buffer, too.
    buffer2 = list()

    sapply(buffer, function(buf) {

      # get the next nodes in this path. NA means no more nodes, i.e no
      # outgoing arc from the last node of the path.
      if (!is.na(buf[length(buf)])) {

        # pick the children of the current node from the adjacency matrix.
        next.one = names(which(amat[buf[length(buf)],] == 1))

        # discard prospective nodes which are already in this path.
        next.one = next.one[!(next.one %in% buf[-1])]

        if (to %in% next.one) {

          # found a path to the "to" node, returning.
          assign("path.found", TRUE, envir = sys.frame(-3))

        }#THEN

        if (length(next.one) >= 1) {

          # if there's someone left, update the path and store it in the
          # temporary buffer ...
          buffer2 = c(buffer2, lapply(next.one, function(x) { c(buf, x) } ))

        }#THEN
        else {

          # ... otherwise terminate the path with a NA so that it will
          # not be updated in the next iterations.
          buffer2 = c(buffer2, list(c(buf, NA)))

        }#ELSE

      }#THEN
      else {

        # drop this path, as there are no new elements.

      }#ELSE

      # update the temporary buffer.
      assign("buffer2", buffer2, envir = sys.frame(-3))

    })

    if (debug) {

      cat("current buffer of discovered paths (depth", i + 1, ") is:\n")
      print(buffer2)

    }#THEN

    # found the path I was looking for, I return TRUE.
    if (path.found) return(TRUE)

    # update the main buffer with the temporary one.
    buffer = buffer2

  }#FOR

  return(FALSE)

}#HAS.PDAG.PATH

# count the loops the arc is part of (in a partially directed graph).
how.many.loops = function(arc, nodes, amat, debug = FALSE) {

  # initialize the loop counter.
  loops = 0

  # if there is no path from arc[2] to arc[1], arc can not be
  # part of a cycle. This fast check (it uses 1-step backtracking,
  # so it may be easily fooled) detects most of the arcs which are
  # not part of any cycle, causing a great performance boost on
  # large graphs.
  if (!has.path(arc[2], arc[1], nodes, amat, exclude.direct = TRUE))
    return(0)

  # initialize the main buffer of the paths.
  buffer = list(arc)

  # loop length(nodes) - 1) times, so that each path has at most
  # 2 + length(nodes) - 1) = length(nodes) + 1 nodes = length(nodes) arcs.
  # no loop can be longer than that unless a node is present more than
  # one time; in that case the loop is dupe of a shorter one (with distinct
  # nodes).
  for (i in 1:(length(nodes) - 1)) {

    # set up a temporary buffer, too.
    buffer2 = list()

    sapply(buffer, function(buf) {

      # get the next nodes in this path. NA means no more nodes, i.e no
      # outgoing arc from the last node of the path.
      if (!is.na(buf[length(buf)])) {

        if (buf[length(buf)] == arc["from"]){

          # if this path is already looping, stop here; the loop
          # counter will be heavily biased if we do otherwise
          # (the same loop, forking over and over and creating
          # lots of fake paths and increasing the loop counter).
          next.one = NA

        }#THEN
        else {

          # pick the children of the current node from the adjacency matrix.
          next.one = names(which(amat[buf[length(buf)],] == 1))

          # discard prospective nodes which are already in the path I'm
          # building; I'm either backtracking, hitting some cycle which
          # does not include the arc I'm interested in at the moment or
          # trying to cross an undirected arc in *both* directons
          # simultaneuosly (ypes!). The obvious exception to the rule:
          # the first element of the buffer, the one which I'm waiting for.
          next.one = next.one[!(next.one %in% buf[-1]) & (next.one != buf[length(buf)-1])]

        }#ELSE

        if (length(next.one) >= 1) {

          # if there's someone left, update the path and store it in the
          # temporary buffer ...
          buffer2 = c(buffer2, lapply(next.one, function(x) { c(buf, x) } ))

        }#THEN
        else {

          # ... otherwise terminate the path with a NA so that it will
          # not be updated in the next iterations.
          buffer2 = c(buffer2, list(c(buf, NA)))

        }#ELSE

      }#THEN
      else {

        # drop this path, as there are no new elements; adn increase the
        # loop counter if there are cycles.
        if (length(which(buf == arc["from"])) > 1) loops = loops + 1

      }#ELSE

      # update the temporary buffer.
      assign("buffer2", buffer2, envir = sys.frame(-3))
      # update the loop counter.
      assign("loops", loops, envir = sys.frame(-3))

    })

    if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("* path discovery for arc", arc["from"], "->", arc["to"], ".\n")
      cat("current buffer is:\n")
      print(buffer2)
      cat("loops the arc is part of:", loops, "\n")
      cat("----------------------------------------------------------------\n")

    }#THEN

    # update the main buffer with the temporary one.
    buffer = buffer2

  }#FOR

  # update and return the loop counter.
  loops + length(which(sapply(buffer,
    function(path) { length(which(path == arc["from"]))}) > 1))

}#HAS.PDAG.PATH

# convert a set of neighbourhoods in an arc list.
mb2arcs = function(mb, nodes) {

  empty.mb = sapply(mb, function(x) {(length(x$nbr) == 0) || is.null(x$nbr) || identical(x$nbr, "")})
  result = do.call(rbind, lapply(nodes[!empty.mb],
               function(x) { cbind(from = x, to = mb[[x]][['nbr']]) }))

  # return an empty matrix all markov blankets are empty.
  if (is.null(result))
    matrix(character(0), ncol = 2, dimnames = list(c(), c("from", "to")))
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

  sapply(nodes(x), nparams.node, x = x, data = data, real = real)

}#NPARAMS.BACKEND

nparams.node = function(node, x, data, real) {

  .Call("nparams",
        graph = x,
        node = node,
        data = data,
        real = as.integer(real),
        PACKAGE = "bnlearn")

}#NPARAMS.NODE

# backend for neighbourhood detection.
nbr.backend = function(arcs, node) {

  # this includes neighbours with undirected arcs.
  unique(c(arcs[arcs[, "from"] == node, "to"], arcs[arcs[, "to"] == node, "from"]))

}#NBR.BACKEND

