
# build a valid whitelist.
build.whitelist = function(whitelist, nodes, data, algo, criterion) {

  if (is.null(whitelist)) {

    # no whitelist, nothing to do.
    return(NULL)

  }#THEN

  if (is(whitelist, c("matrix", "data.frame"))) {

    if (dim(whitelist)[2] != 2)
      stop("whitelist must have two columns.")

    if (is.data.frame(whitelist))
      whitelist = as.matrix(cbind(as.character(whitelist[, 1]),
        as.character(whitelist[, 2])))

  }#THEN
  else if (is.character(whitelist)) {

    if (length(whitelist) != 2)
      stop("whitelist must have two columns.")

    whitelist = matrix(whitelist, ncol = 2, byrow = TRUE)

  }#THEN
  else {

    stop("whitelist must be a matrix or data.frame with two columns.")

  }#ELSE

  # drop duplicate rows.
  whitelist = unique.arcs(whitelist, nodes, warn = TRUE)
  # add column names for easy reference.
  colnames(whitelist) = c("from", "to")

  # check all the names in the whitelist against the column names of x.
  if (any(unique(as.vector(whitelist)) %!in% nodes))
    stop("unknown node label present in the whitelist.")
  # check that whitelisted arcs do not violate parametric assumptions.
  whitelist = check.arcs.against.assumptions(whitelist, data, criterion)

  if (algo %in% score.based.algorithms) {

    # the whitelist should contain only directed arcs.
    if (any(which.undirected(whitelist, nodes = nodes)))
      stop("score-based algorithms do not support whitelisting both directions of an arc.")

  }#THEN
  else if (algo %in% mim.based.algorithms) {

    # all arcs in the whitelist are treated as undirected, because these
    # algorithms operate in the space of undirected graphs.
    whitelist = unique.arcs(arcs.rbind(whitelist, whitelist,
                  reverse2 = TRUE), nodes)

  }#THEN

  # if the whitelist itself contains cycles, no acyclic graph
  # can be learned.
  if (!is.acyclic(whitelist, nodes = nodes,
         directed = (algo %in% c(constraint.based.algorithms, "aracne"))))
    stop("this whitelist does not allow an acyclic graph.")

  return(whitelist)

}#BUILD.WHITELIST

check.arcs.against.assumptions = function(arcs, data, criterion) {

  if (is.null(criterion))
    return(arcs)

  if (criterion %in% c(available.mixedcg.tests, available.mixedcg.scores)) {

    # arcs cannot point from continuous nodes to discrete nodes.
    if (is.null(arcs)) {

      arcs = list.cg.illegal.arcs(nodes = names(data), variables = data)

    }#THEN
    else {

      arcs = .Call(call_arcs_cg_assumptions,
                   arcs = arcs,
                   nodes = names(data),
                   data = data)

    }#ELSE

  }#THEN

  return(arcs)

}#CHECK.ARCS.AGAINST.ASSUMPTIONS

list.cg.illegal.arcs = function(nodes, variables) {

  .Call(call_cg_banned_arcs,
        nodes = nodes,
        variables = variables)

}#LIST.ILLEGAL.ARCS

# build a valid blacklist.
build.blacklist = function(blacklist, whitelist, nodes, algo) {

  if (!is.null(blacklist)) {

    if (is(blacklist, c("matrix", "data.frame"))) {

      if (dim(blacklist)[2] != 2)
        stop("blacklist must have two columns.")

      if (is.data.frame(blacklist))
        blacklist = as.matrix(cbind(as.character(blacklist[, 1]),
          as.character(blacklist[, 2])))

    }#THEN
    else if (is.character(blacklist)) {

      if (length(blacklist) != 2)
        stop("blacklist must have two columns.")

      blacklist = matrix(blacklist, ncol = 2, byrow = TRUE)

    }#THEN
    else {

      stop("blacklist must be a matrix or data.frame with two columns.")

    }#ELSE

    # check all the names in the blacklist against the column names of x.
    if (any(unique(as.vector(blacklist)) %!in% nodes))
      stop("unknown node label present in the blacklist.")

    if (algo %in% mim.based.algorithms) {

      # all arcs in the blacklist are treated as undirected, because these
      # algorithms operate in the space of undirected graphs.
      blacklist = arcs.rbind(blacklist, blacklist, reverse2 = TRUE)

    }#THEN

    # drop duplicate rows.
    blacklist = unique.arcs(blacklist, nodes)

  }#THEN

  # update blacklist to agree with whitelist.
  # NOTE: whitelist and blacklist relationship is the same as hosts.allow
  # and hosts.deny.
  if (!is.null(whitelist)) {

    # if x -> y is whitelisted but y -> x is not, it is to be blacklisted.
    to.add = apply(whitelist, 1, function(x)
               is.whitelisted(whitelist, x[c(2, 1)]))
    blacklist = arcs.rbind(blacklist, whitelist[!to.add, c(2, 1)])

    # if x -> y is whitelisted, it is to be removed from the blacklist.
    if (!is.null(blacklist)) {

      blacklist = blacklist[!apply(blacklist, 1,
        function(x){ is.whitelisted(whitelist, x) }),]

      # also drop duplicate rows.
      blacklist = unique.arcs(matrix(blacklist, ncol = 2, byrow = FALSE), nodes)

    }#THEN

  }#THEN

  # set the column names.
  if (!is.null(blacklist))
    colnames(blacklist) = c("from", "to")

  return(blacklist)

}#BUILD.BLACKLIST

