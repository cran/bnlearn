
# frontend for formulas that count various types of graphs.
count.graphs = function(type = "all.dags", nodes, ..., debug = FALSE) {

  # check whether gmp is loaded, and try to load if it is not.
  check.and.load.package("gmp")
  # check the label describing which graphs to count.
  check.label(type, choices = available.enumerations, argname = "graph types",
    see = "count.graphs")
  # check the number of nodes.
  if (!is.positive.vector(nodes))
    stop("'nodes' must be positive integers, the number(s) of nodes in the graph.")
  largest.graph.size = max(nodes)
  # check debug.
  check.logical(debug)

  # extract the additional arguments.
  extra.args = check.enumeration.args(type, N = nodes, extra = list(...))

  if (type == "all-dags") {

    # count all directed acyclic graphs with N nodes.
    how.many = count.all.dags(N = largest.graph.size, debug = debug)

  }#THEN
  else if (type == "dags-disregarding-one-arc") {

    # count the number of directed acyclic graphs that are identical but for
    # the state of one possible arc between a set of labelled nodes.
    how.many = count.disregarding.one.arc(N = largest.graph.size, debug = debug)

  }#THEN
  else if (type == "dags-given-ordering") {

    # count the number of directed acyclic graphs that satisfy a given
    # topological ordering.
    how.many = count.by.topological.ordering(N = largest.graph.size, debug = debug)

  }#THEN
  else if (type == "dags-with-k-roots") {

    # count the number of directed acyclic graphs that have k roots.
    how.many = count.by.roots(N = largest.graph.size, k = extra.args$k, debug = debug)

  }#THEN
  else if (type == "dags-with-r-arcs") {

    # count the number of directed acyclic graphs that have r arcs.
    how.many = count.by.arcs(N = largest.graph.size, r = extra.args$r, debug = debug)

  }#THEN

  return(how.many[nodes + 1])

}#COUNT.GRAPHS

# count the number of directed acyclic graphs that have r arcs.
count.by.arcs = function(N, r, debug = FALSE) {

  # helper for the maximum possible number of arcs.
  max.arcs = function(nodes) choose(nodes, 2)

  # nothing to do if the number of arcs is too large.
  if (r > max.arcs(N))
    return(gmp::as.bigz(rep(0, N + 1)))

  # allocate a list to hold all the intermediate results, the a[N - k][n]
  # counts of the lower order directed acyclic graphs with various numbers
  # of roots.
  a = vector(N, mode = "list")

  for (i in seq(from = 1, to = N)) {

    # counting graphs with i nodes.
    if (debug)
      cat("* computing the number of DAGs with", i, "nodes.\n")

    # for each graph size, count up to the required number of arcs or the
    # maximum possible number of arcs, whichever is larger.
    a[[i]] = gmp::as.bigz(rep(0, max(max.arcs(i), r) + 1))
    # there is only one graph with zero arcs, the empty graph.
    a[[i]][1] = gmp::as.bigz(1)
    # nothing to do for graphs with one node.
    if (i == 1)
      next

    # compute all the lower-order counts a[i][j].
    for (j in seq(from = 1, to = min(r, max.arcs(i)))) {

      if (debug)
        cat("  >", i, "nodes and", j, "arcs: ")

      # reset the accumulator for the summation over the graphs.
      cumsum = gmp::as.bigz(0)

      # skip the case in which i == m, counts are always zero because they
      # involve graphs with zero nodes.
      for (m in seq(from = 1, to = i - 1)) {

        sign = ifelse((m - 1) %% 2 == 0, +1, -1)
        combinations = gmp::chooseZ(i, m)

        # reset the accumulator of the nested summation over the subgraphs.
        recursion = gmp::as.bigz(0)

        for (k in seq(from = 0, to = min(j, max.arcs(i - m)))) {

          sub.combinations = gmp::chooseZ(m * (i - m), j - k)
          sub.graphs = a[[i - m]][k + 1]
          recursion = recursion + sub.combinations * sub.graphs

        }#FOR

        cumsum = cumsum + sign * combinations * recursion

      }#FOR

      a[[i]][j + 1] = cumsum

      if (debug)
        cat(as.character(a[[i]][j + 1]), "\n")

    }#FOR

  }#FOR

  # now create an array containing all the counts for r arcs; if a graph is too
  # small to have r arcs, the count is missing and the list of intermediate
  # results should not be accessed at all.
  counts = gmp::as.bigz(rep(0, N + 1))
  for (i in seq(from = 1, to = N))
    if (length(a[[i]]) >= r + 1)
      counts[i + 1] = a[[i]][r + 1]

  return(counts)

}#COUNT.BY.ARCS

# count the number of directed acyclic graphs that have k roots.
count.by.roots = function(N, k, debug = FALSE) {

  # nothing to do if the number of roots is too large.
  if (k > N)
    return(gmp::as.bigz(rep(0, N + 1)))

  # allocate a list to hold all the intermediate results, the a[N - k][n]
  # counts of the lower order directed acyclic graphs with various numbers
  # of roots.
  a = vector(N, mode = "list")

  for (i in seq(from = 1, to = N)) {

    # counting graphs with i nodes.
    if (debug)
      cat("* computing the number of DAGs with", i, "nodes.\n")

    a[[i]] = gmp::as.bigz(rep(0, i))

    # compute all the lower-order counts a[i][j].
    for (j in seq(from = 1, to = max(i, 1))) {

      if (debug)
        cat("  > with", j, "roots: ")

      if (i == j) {

        # if all nodes are roots, the only possible graph is the empty graph.
        a[[i]][j] = gmp::as.bigz(1)

      }#THEN
      else {

        combinations = gmp::chooseZ(i, j)

        for (n in seq(from  = 1, to = i - j)) {

          arcs.to.old.roots = gmp::pow.bigz(gmp::pow.bigz(2, j) - 1, n)
          arcs.to.old.nonroots = gmp::pow.bigz(2, j * (i - n - j))
          recurse = a[[i - j]][n]

          a[[i]][j] = a[[i]][j] + arcs.to.old.roots * arcs.to.old.nonroots * recurse

        }#FOR

        a[[i]][j] = a[[i]][j] * combinations

      }#ELSE

      if (debug)
        cat(as.character(a[[i]][j]), "\n")

    }#FOR

  }#FOR

  # now create an array containing all the counts for k roots; the counts are
  # automatically zero if the number of roots is larger than the numnber of
  # nodes.
  counts = gmp::as.bigz(rep(0, N + 1))
  for (i in seq(from = k, to = N))
    counts[i + 1] = a[[i]][k]

  return(counts)

}#COUNT.BY.ROOTS

# count all directed acyclic graphs with N nodes.
count.all.dags = function(N, debug = FALSE) {

  # allocate the return value, using extended precision.
  a = gmp::as.bigz(rep(NA, N + 1))
  # initialize the sequence.
  a[1] = 1

  for (i in seq(from = 1, to = N)) {

    # computing a_i, storing in a[i + 1] because a_0 is in a[1].
    if (debug)
      cat("* computing the number of DAGs with", i, "nodes: ")

    # reset the summation accumulator.
    cumsum = gmp::as.bigz(0)

    # summation in the recusrive formula.
    for (k in i:1) {

      sign = ifelse((k - 1) %% 2 == 0, +1, -1)
      combinations = gmp::chooseZ(i, k)
      possible.arcs = gmp::pow.bigz(2, k * (i - k))

      cumsum = cumsum + sign * combinations * possible.arcs * a[i - k + 1]

    }#FOR

    # save the results in the return value.
    a[i + 1] = cumsum

    if (debug)
      cat(as.character(a[i + 1]), "\n")

  }#FOR

  return(a)

}#COUNT.ALL.DAGS

# count the number of directed acyclic graphs that are identical but for
# the state of one possible arc between a set of labelled nodes.
count.disregarding.one.arc = function(N, debug = FALSE) {

  # first, we need to compute all the a[i + 1].
  a = count.all.dags(N - 1, debug = debug)
  # then, allocate the return value and initialize the sequence.
  b = gmp::as.bigz(rep(NA, N + 1))
  b[1:2] = 1

  for (i in seq(from = 2, to = N)) {

    # computing b_i, storing in b[i + 1] so that b_0 is in b[1].
    if (debug)
      cat("* computing the number of DAGs with", i, "nodes, modulo one arc: ")

    # reset the summation accumulator.
    cumsum = gmp::as.bigz(0)

    # summation in the recusrive formula.
    for (k in i:1) {

      sign = ifelse((k - 1) %% 2 == 0, +1, -1)
      combinations = gmp::chooseZ(i - 1, k - 1)
      possible.arcs = gmp::pow.bigz(2, k * (i - k))
      reduced.combinations = gmp::chooseZ(i - 2, k)

      cumsum = cumsum + sign * possible.arcs *
        (reduced.combinations * b[i - k + 1] + combinations * a[i - k + 1])

    }#FOR

    b[i + 1] = cumsum

    if (debug)
      cat(as.character(b[i + 1]), "\n")

  }#FOR

  return(b)

}#COUNT.DISREGARDING.ONE.ARC

# count the number of directed acyclic graphs that satisfy a given
# topological ordering.
count.by.topological.ordering = function(N, debug = FALSE) {

  # allocate the return value and initialize the sequence.
  a = gmp::as.bigz(rep(NA, N + 1))
  a[1:2] = 1

  for (i in seq(from = 2, to = N)) {

    # computing a_i, storing in a[i + 1] so that a_0 is in b[1].
    if (debug)
      cat("* computing the number of DAGs with", i, "nodes, one arc fixed.\n")

    a[i + 1] = gmp::pow.bigz(2, gmp::chooseZ(i, 2))

  }#FOR

  return(a)

}#COUNT.BY.TOPOLOGICAL.ORDERING
