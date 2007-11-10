
# get all the subsets of a given size, even if either the initial set
# or the subset are empty (i.e. of size zero).
subsets = function (n, r, v = 1:n, set = TRUE, repeats.allowed = FALSE) {

  # allow empty subsets (i.e. subsets of empty sets).
  if ((n == 0) || (r == 0)) return(matrix(c(""),1,1))

  if (mode(n) != "numeric" || length(n) != 1 || n < 1 || (n%%1) != 0)
    stop("bad value of n")

  if (mode(r) != "numeric" || length(r) != 1 || r < 1 || (r%%1) != 0)
    stop("bad value of r")

  if (!is.atomic(v) || length(v) < n)
    stop("v is either non-atomic or too short")

  if ((r > n) & repeats.allowed == FALSE)
    stop("r > n and repeats.allowed=FALSE")

  if (set) {

    v <- unique(sort(v))
    if (length(v) < n)
      stop("too few different elements")

  }#THEN

  v0 <- vector(mode(v), 0)
  if (repeats.allowed) {

    sub <- function(n, r, v) {

        if (r == 0)
            v0
        else if (r == 1)
            matrix(v, n, 1)
        else if (n == 1)
            matrix(v, 1, r)
        else rbind(cbind(v[1], Recall(n, r - 1, v)), Recall(n - 1, r, v[-1]))

    }#SUB

  }#THEN
  else {

    sub <- function(n, r, v) {

        if (r == 0)
            v0
        else if (r == 1)
            matrix(v, n, 1)
        else if (r == n)
            matrix(v, 1, n)
        else rbind(cbind(v[1], Recall(n - 1, r - 1, v[-1])), Recall(n - 1, r, v[-1]))

    }#SUB

  }#ELSE

  sub(n, r, v[1:n])

}#SUBSETS

# check a matrix for symmetry.
is.symmetric = function(m) {

  # kill all the names; identical may return false otherwise.
  colnames(m) = rownames(m) = NULL

  identical(m, t(m))

}#IS.SYMMETRIC

# check which rows of a data frame or matrix are identical to an array.
is.row.equal = function(data, array) {

  apply(data, 1, function(x){ all(x == array) })

}#IS.ROW.EQUAL

# return the array whose size is smaller.
smaller = function(a, b) {

  if (length(a) < length(b))
    a
  else
    b

}#SMALLER

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
is.undirected = function(arcs) {

  apply(arcs, 1, function(arc) { is.listed(arcs, arc, both = TRUE) })

}#WHICH.UNDIRECTED

# set a direction for an (undirected) arc
set.arc.direction = function(from, to, arcs, debug = FALSE) {

  # the arc is there, undirected
  if (is.listed(arcs, c(to, from), both = TRUE)) {

    if (debug)
      cat("  > the arc", from, "-", to, "is undirected,",
            "changing to", from, "->", to, ".\n")

    arcs[!is.row.equal(arcs, c(to, from)),]

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
         is.row.equal(arcs, dropped[c(2,1)])),]

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

# is there a path from a to b (besides a -> b, if that's the case)?
has.path = function(arcs, nodes, a, b) {

   # use the adjacency matrix as a non-stochastic transition matrix
   # in a markov chain.
   amat = arcs2amat(arcs, nodes, stochastic = TRUE)
   start = as.integer(nodes == a)

   for (i in 1:(length(nodes)+1)) {

     # the multiplication must be done before the check; otherwise
     # the return value is always true.
     amat = amat %*% amat

     if ((t(start) %*% amat)[, b] > 0) return(TRUE)

   }#FOR

   FALSE

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

# check whether a graph is acyclic.
is.acyclic = function(x) {

  all(loop.counter(x$arcs, names(x$nodes))[,3] == 0)

}#IS.ACYCLIC

# compute the data / cells ratio.
obs.per.cell = function(x, y, z = NULL, data) {

  if (is.null(z) || (length(z) == 0)) {

    nrow(data) / (nlevels(data[,x]) * nlevels(data[,y]))

  }#THEN
  else if (is.character(z)) {

    if (length(z) == 1)
      nrow(data) / (nlevels(data[,x]) * nlevels(data[,y]) * nlevels(data[,z]))
    else if (length(z) > 1)
      nrow(data) / (nlevels(data[,x]) * nlevels(data[,y]) *
        prod(sapply(z, function(col) { nlevels(data[, col]) } )))

  }#THEN
  else if (is.factor(z)) {

    nrow(data) / (nlevels(data[,x]) * nlevels(data[,y]) * nlevels(z))

  }#ELSE

}#OBS.PER.CELL


# convert a set of markov blankets to a (sort of) adjacency matrix.
mb2amat = function(mb, stochastic = FALSE, redux = FALSE) {

  if (redux)
    amat = sapply(names(mb), function(x) {as.numeric(names(mb) %in% mb[[x]])})
  else
    amat = sapply(names(mb), function(x) {as.numeric(names(mb) %in% mb[[x]]$mb)})

  if (stochastic)
    amat = t(apply(amat, 1, function(r) { if(sum(r) != 0) r / sum(r) else r}))

  amat

}#MB2AMAT

# convert a set of neighbourhoods to a (real) adjacency matrix.
nbr2amat = function(mb, stochastic = FALSE) {

  amat = sapply(names(mb), function(x) {as.numeric(names(mb) %in% mb[[x]]$nbr)})

  if (stochastic)
    amat = t(apply(amat, 1, function(r) { if(sum(r) != 0) r / sum(r) else r}))

  amat

}#NBR2AMAT

# convert a set of arcs to a (real) adjacency matrix.
arcs2amat = function(arcs, nodes, stochastic = FALSE) {

  amat = matrix(rep(0, length(nodes)^2), ncol = length(nodes))
  colnames(amat) = rownames(amat) = nodes

  if (nrow(arcs) > 0)
    for (i in 1:nrow(arcs)) {

      amat[arcs[i, "from"], arcs[i, "to"]] = 1

    }#FOR

  if (stochastic)
    amat = t(apply(amat, 1, function(r) { if(sum(r) != 0) r / sum(r) else r}))

  amat

}#ARCS2AMAT

# convert an adjacency matrix back to a set of arcs.
amat2arcs = function(a, nodes, debug = FALSE) {

  if (is.null(colnames(a)))
    colnames(a) = nodes
  if (is.null(rownames(a)))
    rownames(a) = nodes

  arcs = do.call(rbind,
    sapply(rownames(a), function(node) {

      tos = colnames(a)[as.logical(a[node, ])]

      if (debug) {

        cat("* preocessing node", node, "\n")
        cat("  > adding arcs from", node, "to its parents: '", tos, "'.\n")

      }#THEN

      if (length(tos) > 0)
        cbind(rep(node, length(tos)), tos)

    })
  )

  colnames(arcs) = c("from", "to")

  arcs

}

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

# check the cluster is running.
isClusterRunning = function(cl) {

  tryCatch(any(clusterEvalQ(cl, TRUE)),
    error = function(err) { FALSE })

}#ISCLUSTERRUNNING

# print an underlined label in the plot.
underlined <- function(x, y, label, col){

  text(x, y, label, col = col, font = 2)
  sw <- strwidth(label)
  sh <- strheight(label)
  lines(x + c(-sw/2, sw/2), rep(y - 1.5*sh/2, 2), col = col)

}#UNDERLINED

# will the bayesian network be a discrete one?
is.data.discrete = function(data) {

  all(sapply(data, class) == "factor")

}#IS.DATA.DISCRETE

# will the bayesian network be a continuous one?
is.data.continuous = function(data) {

  all(sapply(data, class) == "numeric")

}#IS.DATA.CONTINUOUS

# there are missing data?
missing.data = function(data) {

  any(mapply(function(x) {is.na(x) || is.nan(x) || is.null(x)}, data))

}#MISSING.DATA

