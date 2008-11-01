# check which rows of a data frame or matrix are identical to an array.
is.row.equal = function(data, array) {

  .Call("is_row_equal",
        data = as.character(data),
        array = as.character(array),
        PACKAGE = "bnlearn")

}#IS.ROW.EQUAL

# check whether an arc is present in a matrix or data frame (with 2 columns).
is.listed = function(list, arc, either = FALSE, both = FALSE) {

  if (both && either)
    stop("conflicting options both and either.")

  if (is.null(list))
    FALSE
  else if (is.null(arc))
    stop(" a valid arc must be specified.")
  else if (both)
    any(is.row.equal(list, arc)) && any(is.row.equal(list, arc[c(2,1)]))
  else if (either)
    any(is.row.equal(list, arc)) || any(is.row.equal(list, arc[c(2,1)]))
  else
    any(is.row.equal(list, arc))

}#IS.LISTED

is.whitelisted = is.listed
is.blacklisted = is.listed

# which arcs are undirected?
which.undirected = function(arcs) {

  .Call("which_undirected",
        arcs = factor(arcs),
        PACKAGE = "bnlearn")

}#WHICH.UNDIRECTED

# is this arc undirected?
is.undirected = function(arc, arcs) {

  # if it is there must be its reverse in the arc set.
  is.listed(arcs, arc[c(2,1)])

}#IS.UNDIRECTED

#is this arc directed?
is.directed = function(arc, arcs) {

  !is.undirected(arc, arcs)

}#IS.DIRECTED

# set a direction for an (undirected) arc
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

    rbind(arcs[!is.row.equal(arcs, c(to, from)),], c(from, to))

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

    rbind(arcs, c(from, to))

  }#ELSE

}#SET.ARC.DIRECTION

# drop an arc.
drop.arc.backend = function(arcs, dropped, debug = FALSE) {

  if (debug)
    cat("  > dropping any arc between", dropped[1], "and", dropped[2], ".\n")

  # drop the arc.
  arcs[!(is.row.equal(arcs, dropped) |
         is.row.equal(arcs, dropped[c(2,1)])), , drop = FALSE]

}#DROP.ARC.BACKEND

# reverse the direction of an arc.
reverse.arc.backend = function(from, to, arcs, debug = FALSE) {

  # the arc is there, undirected
  if (is.listed(arcs, c(to, from), both = TRUE))
    stop("an undirected arc cannot be reversed.")
  # the arc is there, but reversed.
  else if (is.listed(arcs, c(to, from))) {

    if (debug)
      cat("  > the arc", to, "->", from, "is present, reversing.\n")

    rbind(arcs[!is.row.equal(arcs, c(to, from)),], c(from, to))

  }#THEN
  # the arc is already there.
  else if (is.listed(arcs, c(from, to))) {

    if (debug)
      cat("  > the arc", from, "->", to, "is present, reversing.\n")

    rbind(arcs[!is.row.equal(arcs, c(from, to)),], c(to, from))

  }#THEN
  # the arc is not present.
  else
    stop("no arc to be reversed between ", from, " and ", to, ".\n")

}#REVERSE.ARC.BACKEND

# which arcs are {white,black}listed?
which.listed = function(arcs, list) {

  apply(arcs, 1, function(arc) { is.listed(list, arc) })

}#WHICH.LISTED

which.whitelisted = which.listed
which.blacklisted = which.listed

