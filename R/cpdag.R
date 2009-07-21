
# reconstruct the equivalence class of a network.
cpdag.backend = function(amat, nodes, debug = TRUE) {

  # convert the adjacency matrix to an object of class bn.
  bn = empty.graph(nodes)
  amat(bn) = amat

  # separate direct ted and undirected arcs.
  uarcs = undirected.arcs(bn)
  darcs = directed.arcs(bn)

  # extract compelled arcs (i.e. those which are part of a v-structure).
  tmp = table(darcs[, "to"])
  vs.child = names(tmp[which(tmp > 1)])
  compelled.arcs = darcs[darcs[, "to"] %in% vs.child , , drop = FALSE]

  other.arcs = bn$arcs[!which.listed(bn$arcs, compelled.arcs) , , drop = FALSE]

  res = empty.graph(nodes)

  # compelled arcs are preserved in the equivalence class, apply.
  arcs(res) = compelled.arcs

  if (debug) {

    cat("all arcs:\n")
    print(bn$arcs)
    cat("other arcs:\n")
    print(other.arcs)
    cat("compelled arcs:\n")
    print(compelled.arcs)

  }#THEN

  to.drop = numeric(0)

  repeat {

    if (nrow(other.arcs) == 0) break

    if (debug) cat("nrow(other.arcs) = ", nrow(other.arcs), "\n")

    for (i in 1:nrow(other.arcs)) {

      if (debug) cat("* considering arc", other.arcs[i, 1], "->", other.arcs[i, 2], ".\n")

      if (is.listed(res$arcs, other.arcs[i, , drop = FALSE], either = TRUE)) {

        if (debug) cat("  > dropping arc", other.arcs[i, 1], "->", other.arcs[i, 2], ".\n")

        # this arc has already been set in the network, dropping.
        to.drop = c(to.drop, i)
        next

      }#THEN

      if (path(res, other.arcs[i, "from"], other.arcs[i, "to"])) {

        if (debug) cat("  > path from", other.arcs[i, 1], "to", other.arcs[i, 2], "setting", other.arcs[i, 1], "->", other.arcs[i, 2], ".\n")

        # if there is a directed path from A to B, drop A -B in favour of A -> B.
        res = set.arc(res, other.arcs[i, "from"], other.arcs[i, "to"])

        to.drop = c(to.drop, i)
        next

      }#THEN

      for (node in nodes[!(nodes %in% other.arcs[i, ])]) {

        # A = other.arcs[i, 1], B = node, C = other.arcs[i, 2]
        # if A is not adjacent to C, but A - B and C -> B, set A -> B.

        if (!is.listed(bn$arcs, c(other.arcs[i, 1], node), either = TRUE) &&
             is.listed(res$arcs, c(node, other.arcs[i, 2])) &&
            !is.listed(res$arcs, other.arcs[i, ], either = TRUE) ) {

          if (debug) {

            cat("  > killing spurious v-structure in", node, "->", other.arcs[i, "to"], "-", other.arcs[i, "from"] ,".\n")
            cat("  > setting", other.arcs[i, "to"], "->", other.arcs[i, "from"], ".\n")

          }#THEN

          res = set.arc(res, other.arcs[i, "to"], other.arcs[i, "from"])
          to.drop = c(to.drop, i)
          next

        }#THEN

        if (!is.listed(bn$arcs, c(other.arcs[i, 2], node), either = TRUE) &&
             is.listed(res$arcs, c(node, other.arcs[i, 1])) &&
            !is.listed(res$arcs, other.arcs[i, c(2,1)], either = TRUE) ) {

          if (debug) {

            cat("  > killing spurious v-structure in", other.arcs[i, "from"], "-", other.arcs[i, "to"], "<-", node, ".\n")
            cat("  > setting", other.arcs[i, "from"], "->", other.arcs[i, "to"], ".\n")

          }#THEN

          res = set.arc(res, other.arcs[i, "from"], other.arcs[i, "to"])
          to.drop = c(to.drop, i)
          next

        }#THEN

      }#FOR

    }#FOR

    if (length(to.drop) >= 1)
      other.arcs = other.arcs[-to.drop, , drop = FALSE]
    else if (length(to.drop) == 0)
      break

    to.drop = numeric(0)

  }#REPEAT

  res$arcs = matrix(c(res$arcs[, 1], other.arcs[, 1], other.arcs[, 2], res$arcs[, 2], 
              other.arcs[, 2], other.arcs[, 1]), ncol = 2, dimnames = list(NULL, c("from", "to")))
  res$arcs = unique(res$arcs)

  return(res)

}#CPDAG

