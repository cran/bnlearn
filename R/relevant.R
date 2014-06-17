
# Pena's irrelevant nodes feature selection.
pena.backend  = function(target, context, data, test, alpha, B, debug = TRUE) {

  to.test = relevant = target
  nodes = setdiff(names(data), context)
  next.to.test = to.add = character(0)
  changed = FALSE

  reset.test.counter()

  if (debug) {

    cat("* target nodes are '", target, "'\n")
    cat("* context nodes are '", context, "'\n")

  }#THEN

  repeat {

    for (node in to.test) {

      if (debug) {

        if (is.null(context))
          cat("* testing for marginal association with", node, ".\n")
        else
          cat("* testing for conditional association with", node, ".\n")

      }#THEN

      # compute the marginal associations.
      association = indep.test(setdiff(nodes, relevant), y = node, sx = context,
                      data = data, test = test, alpha = alpha, B = B)

     if (debug) {

       sapply(names(association),
          function(x) {  cat("    >", x, "has p-value", association[x], ".\n")})

      }#THEN

      # select which nodes to add according to the alpha threshold.
      to.add = names(which(association < alpha))
      # check whether there are new relevant nodes.
      changed = changed || length(to.add > 0)
      # update the set of relevant nodes.
      relevant = union(relevant, to.add)
      # update the next set of variables to test.
      next.to.test = union(next.to.test, to.add)

      if (debug && (length(to.add) > 0)) {

        cat("  > adding nodes '", to.add, "'\n")
        cat("  @ relevant nodes are '", relevant, "'\n")

      }#THEN

    }#FOR

    if (!changed) {

      # no new node has been added, we are done.
      break

    }#THEN
    else {

      # reset the changed variable.
      changed = FALSE
      # update the set of nodes to test.
      to.test = next.to.test
      # reset the next set of variables to test.
      next.to.test = character(0)

    }#ELSE

  }#REPEAT

  # remove the target variables and return.
  return(setdiff(relevant, target))

}#PENA.BACKEND

