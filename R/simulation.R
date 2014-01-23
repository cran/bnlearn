# do a partial ordering of the nodes of a graph.
schedule = function(x, start = NULL, reverse = FALSE, debug = FALSE) {

  if (is.null(start))
    nodes = root.leaf.nodes(x, leaf = reverse)
  else
    nodes = start

  to.do = .Call("schedule",
                bn = x,
                root.nodes = nodes,
                reverse = reverse,
                debug = debug)

  if (is.null(start))
    return(names(sort(to.do)))
  else
    return(names(sort(to.do[to.do > 0])))

}#SCHEDULE

# use the Logic Sampling (LS) algorithm as described in "Bayesian Artificial
# Intelligence", Korb & Nicholson, chap 3.6.1.
rbn.backend = function(x, n, data, fix = TRUE, debug = FALSE) {

  # fit the bayesian network if needed.
  if (is(x, "bn"))
    fitted = bn.fit.backend(x, data, debug = FALSE)
  else
    fitted = x

  .Call("rbn_master",
        fitted = fitted,
        n = as.integer(n),
        fix = fix,
        debug = debug)

}#RBN.BACKEND

