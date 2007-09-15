
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
which.undirected = function(arcs) {

  apply(arcs, 1, function(arc) { is.listed(arcs, arc, both = TRUE) })

}#IS.UNDIRECTED

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

# get the parents of a node
parents.backend = function(arcs, node, undirected = FALSE) {

  if (!undirected)
    arcs[(arcs[, "to"] == node) & !which.undirected(arcs), "from"]
  else
    arcs[(arcs[, "to"] == node), "from"]

}#PARENTS.BACKEND

loglik.node = function(node, x, data, debug) {

  if (debug) cat("* processing node", node, ".\n")

  # this node is a root node.
  if (length(parents(x, node)) == 0) {

    ndata = nrow(data)
    tab = table(data[, node])
    node.loglik = log(dmultinom(tab, ndata, tab / ndata))

  }#THEN
  # this node has at least one parent.
  else {

    config = factor(apply(as.data.frame(data[, parents(x, node)]), 
      1, paste, sep = "", collapse = ":"))

    tab = table(data.frame(node = data[, node], config))
    node.loglik = log(sum(apply(tab, 2, 
      function(t){ dmultinom(t, sum(t), t / sum(t)) })))

  }#ELSE

  if (debug) { 

    cat("  > node contribution is", node.loglik, ".\n")

  }#THEN

  node.loglik

}#LOGLIK.NODE

# check the cluster is running.
isClusterRunning = function(cl) {

  tryCatch(any(clusterEvalQ(cl, TRUE)),
    error = function(err) { FALSE })

}#ISCLUSTERRUNNING
