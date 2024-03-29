# check which rows of a data frame or matrix are identical to an array.
is.row.equal = function(data, array) {

  .Call(call_is_row_equal,
        data = as.character(data),
        array = as.character(array))

}#IS.ROW.EQUAL

# check whether an arc is present in a matrix or data frame (with 2 columns).
is.listed = function(set, arc, either = FALSE, both = FALSE) {

  .Call(call_is_listed,
        arc = as.character(arc),
        set = as.character(set),
        either = either,
        both = both,
        debug = FALSE)

}#IS.LISTED

is.whitelisted = is.listed
is.blacklisted = is.listed

# which arcs are undirected?
which.undirected = function(arcs, nodes = NULL) {

  .Call(call_which_undirected,
        arcs = arcs,
        nodes = nodes)

}#WHICH.UNDIRECTED

# which arcs are directed?
which.directed = function(arcs, nodes = NULL) {

  !which.undirected(arcs, nodes = nodes)

}#WHICH.DIRECTED

# is this arc undirected?
is.undirected = function(arc, arcs) {

  # both A -> B  and B -> A must be present in the arc set.
  is.listed(arcs, arc, both = TRUE)

}#IS.UNDIRECTED

# is this arc directed?
is.directed = function(arc, arcs) {

  !is.undirected(arc, arcs)

}#IS.DIRECTED

# set the direction of an (un)directed arc.
set.arc.direction = function(from, to, arcs, debug = FALSE) {

  # the arc is there, undirected
  if (is.listed(arcs, c(to, from), both = TRUE)) {

    if (debug)
      cat("  > the arc", from, "-", to, "is undirected,",
            "changing to", from, "->", to, ".\n")

    arcs[!is.row.equal(arcs, c(to, from)), , drop = FALSE]

  }#THEN
  # the arc is there, but the direction is wrong.
  else if (is.listed(arcs, c(to, from))) {

    if (debug)
      cat("  > the arc", to, "->", from, "is present, reversing.\n")

    arcs.rbind(arcs[!is.row.equal(arcs, c(to, from)),], c(from, to))

  }#THEN
  # the arc is already there.
  else if (is.listed(arcs, c(from, to))) {

    if (debug)
      cat("  > the arc", from, "->", to, "is already present, nothing to do.\n")

    arcs

  }#THEN
  # the arc is not present.
  else {

    if (debug)
      cat("  > the arc", from, "->", to, "is not present, adding.\n")

    arcs.rbind(arcs, c(from, to))

  }#ELSE

}#SET.ARC.DIRECTION

# set an undirected edge.
set.edge.backend = function(from, to, arcs, debug = FALSE) {

  # the arc is there, undirected
  if (is.listed(arcs, c(to, from), both = TRUE)) {

    if (debug)
      cat("  > the arc", from, "-", to, "is undirected, nothing to do.\n")

    arcs

  }#THEN
  # the arc is there, but the direction is wrong.
  else if (is.listed(arcs, c(to, from))) {

    if (debug)
      cat("  > the arc", to, "->", from, "is present, removing direction.\n")

    arcs.rbind(arcs, c(from, to))

  }#THEN
  # the arc is already there.
  else if (is.listed(arcs, c(from, to))) {

    if (debug)
      cat("  > the arc", from, "->", to, "is present, removing direction.\n")

    arcs.rbind(arcs, c(to, from))

  }#THEN
  # the arc is not present.
  else {

    if (debug)
      cat("  > the arc", from, "-", to, "is not present, adding.\n")

    arcs.rbind(arcs, c(from, to, to, from))

  }#ELSE

}#SET.EDGE

# drop an arc.
drop.arc.backend = function(arcs, dropped, debug = FALSE) {

  # the arc is there, undirected.
  if (debug)
    cat("  > dropping any arc between", dropped[1], "and", dropped[2], ".\n")

  arcs[!(is.row.equal(arcs, dropped) |
         is.row.equal(arcs, dropped[c(2, 1)])), , drop = FALSE]

}#DROP.ARC.BACKEND

# drop an edge.
drop.edge.backend = function(arcs, dropped, debug = FALSE) {

  # the arc is there, undirected.
  if (is.listed(arcs, dropped, both = TRUE)) {

    if (debug)
      cat("  > dropping any arc between", dropped[1], "and", dropped[2], ".\n")

    arcs[!(is.row.equal(arcs, dropped) |
           is.row.equal(arcs, dropped[c(2, 1)])), , drop = FALSE]

  }#THEN
  else {

    if (debug)
      cat("  > the arc", dropped[1], "-", dropped[2], "is not present, nothing to do.\n")

    arcs

  }#ELSE

}#DROP.EDGE.BACKEND

# reverse the direction of an arc.
reverse.arc.backend = function(from, to, arcs, debug = FALSE) {

  # the arc is there, undirected.
  if (is.listed(arcs, c(to, from), both = TRUE))
    stop("an undirected arc cannot be reversed.")
  # the arc is there, but reversed.
  else if (is.listed(arcs, c(to, from))) {

    if (debug)
      cat("  > the arc", to, "->", from, "is present, reversing.\n")

    arcs.rbind(arcs[!is.row.equal(arcs, c(to, from)),], c(from, to))

  }#THEN
  # the arc is already there.
  else if (is.listed(arcs, c(from, to))) {

    if (debug)
      cat("  > the arc", from, "->", to, "is present, reversing.\n")

    arcs.rbind(arcs[!is.row.equal(arcs, c(from, to)),], c(to, from))

  }#THEN
  # the arc is not present.
  else
    stop("no arc to be reversed between ", from, " and ", to, ".\n")

}#REVERSE.ARC.BACKEND

# which arcs are {white,black}listed?
which.listed = function(arcs, list, either = FALSE, both = FALSE) {

  apply(arcs, 1, function(arc) { is.listed(list, arc, either, both) })

}#WHICH.LISTED

# convert a set of neighbourhoods to an arc set.
nbr2arcs = function(nbr) {

  .Call(call_nbr2arcs,
        nbr = nbr)

}#NBR2ARCS

# remove duplicate arcs and re-orient them according to the node ordering
# specified by the labels.
arcs.unique = function(arcs, nodes, warn = FALSE) {

  .Call(call_unique_arcs,
        arcs = arcs,
        nodes = nodes,
        warn = warn)

}#ARCS.UNIQUE

# rbind-like function for arc sets.
arcs.rbind = function(matrix1, matrix2, reverse2 = FALSE) {

  .Call(call_arcs_rbind,
        matrix1 = matrix1,
        matrix2 = matrix2,
        reverse2 = reverse2)

}#ARCS.RBIND

# return the arcs from an object of class bn.fit.
fit2arcs = function(x) {

  .Call(call_fit2arcs,
        x = x)

}#FIT2ARCS

# return the size of the arc set from an object of class bn or bn.fit.
narcs.backend = function(x) {

  .Call(call_num_arcs,
        x = x)

}#NARCS.BACKEND
