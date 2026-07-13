
# reconstruct the equivalence class of a network.
cpdag.backend = function(x, moral = FALSE, fix = FALSE, wlbl = TRUE,
    debug = FALSE) {

  nodes = names(x$nodes)

  # reset wlbl if x contains no whitelist and no blacklist.
  if ((is.null(x$learning$whitelist)) && (is.null(x$learning$blacklist)))
    wlbl = FALSE

  amat = .Call(call_cpdag,
               arcs = x$arcs,
               nodes = nodes,
               moral = moral,
               fix = fix,
               wlbl = wlbl,
               whitelist = x$learning$whitelist,
               blacklist = x$learning$blacklist,
               illegal = x$learning$illegal,
               debug = debug)

  # update the arcs of the network.
  x$arcs = amat2arcs(amat, nodes)
  # update the network structure.
  x$nodes = cache.structure(nodes, amat = amat)

  return(x)

}#CPDAG.BACKEND

# backend to get a DAG out of a CPDAG (still in the same equivalence class).
cpdag.extension = function(x, debug = FALSE) {

  nodes = names(x$nodes)

  # update the arcs of the network with a consistent extension.
  x$arcs = .Call(call_pdag_extension,
                 arcs = x$arcs,
                 nodes = nodes,
                 debug = debug)
  # update the network structure.
  x$nodes = cache.structure(nodes, arcs = x$arcs)

  return(x)

}#CPDAG.EXTENSION

# produce all possible extensions of a CPDAG.
cextend.all.backend = function(x, debug = FALSE) {

  # if all arcs are directed, there is nothing to extend.
  nodes = names(x$nodes)
  undirected = which.undirected(x$arcs, nodes)

  if (!any(undirected))
    return(x)

  # drop all directed arcs.
  ug = empty.graph.backend(nodes)
  ug$arcs = x$arcs[undirected, , drop = FALSE]
  ug$nodes = cache.structure(nodes, arcs = ug$arcs)

  # drop all isolated nodes, which are irrelevant to the algorithm.
  isolated = (sapply(names(ug$nodes), .degree, x = ug) == 0)
  ug = subgraph.backend(ug, names(which(!isolated)))

  # reset the recursion depth limit to allow the expansion.
  recursion.depth.limit = options("expressions")
  if (recursion.depth.limit < length(ug$nodes)) {

    # the largest valid value is 500000.
    options("expressions" = min(length(ug$nodes) + 50, 500000))
    on.exit(options("expressions" = recursion.depth.limit))

  }#THEN

  # produce all possible combinations of directions of the undirected arcs.
  extensions = uccg.all.extensions(ug, debug = debug)

  # put back the directed arcs.
  extensions = lapply(extensions, function(ext) {

    dag = empty.graph.backend(nodes)
    dag$arcs = rbind(ext, x$arcs[!undirected, , drop = FALSE])
    dag$nodes = cache.structure(nodes, arcs = dag$arcs)

    return(dag)

  })

  return(extensions)

}#CEXTEND.ALL.BACKEND

# produce all consistent extensions of a PDAG or MPDAG, following the bucket
# generalisation of the MCS enumeration algorithm in Wienobst, Luttermann,
# Bannach and Liskiewicz (2023), "Efficient Enumeration of Markov Equivalent
# DAGs" (AAAI 37:12313-12320), Section 4 and Algorithm 7. This backend is
# deliberately kept separate from cextend.all.backend(), which handles the
# simpler CPDAG case where the directed arcs are all compelled.
#
# work in progress: steps C2 (maximal orientation) and C3 (bucket extraction)
# are implemented below; the S+-restricted AMO enumeration (C4) and the
# re-attachment of the between-bucket arcs (C5) are not done yet.
cextend.all.pdag.backend = function(x, debug = FALSE) {

  nodes = names(x$nodes)

  # C6: reject non-extendable inputs up front. A PDAG has a consistent extension
  # if and only if a Dor-Tarsi orientation directs it completely; equivalently
  # (Wienobst et al. 2023, Fact 3) if and only if every bucket skeleton is
  # chordal. Testing it here, before the maximal orientation, also rules out
  # PDAGs whose Meek closure would otherwise force a directed cycle.
  if (!is.completely.directed(cpdag.extension(x))) {

    if (debug)
      cat("* 'x' has no consistent extension.\n")

    stop("'x' is not an extendable PDAG: it has no consistent extension.")

  }#THEN

  # C2: maximally orient the PDAG into its MPDAG with Meek's rules, preserving
  # the directed arcs given as background knowledge (Algorithm 7, line 4).
  amat = arcs2amat(x$arcs, nodes)
  amat = meek.rules.closure(amat, nodes, debug = debug)

  if (debug) {

    cat("* maximally-oriented PDAG (MPDAG):\n")
    print(amat2arcs(amat, nodes))

  }#THEN

  # C3: split the MPDAG into buckets and classify the directed arcs into those
  # that live inside a bucket and those that connect two different buckets
  # (Definition 1; Algorithm 7, lines 6-7).
  buckets = pdag.buckets(amat, nodes, debug = debug)

  # C4: enumerate the AMOs of the buckets with the S+ restriction, so that every
  # ordering conforms to the within-bucket directed arcs (Algorithm 6).
  amos = pdag.buckets.amos(buckets, nodes, debug = debug)

  if (debug)
    cat("* enumerated", length(amos), "bucket AMO(s).\n")

  # C5: re-attach the between-bucket directed arcs to every AMO and wrap each
  # extension into a bn object (Algorithm 7, line 7).
  extensions = lapply(amos, function(amo) {

    dag = empty.graph.backend(nodes)
    dag$arcs = rbind(amo, buckets$between.arcs)
    dag$nodes = cache.structure(nodes, arcs = dag$arcs)

    return(dag)

  })

  if (debug)
    cat("* produced", length(extensions), "consistent extension(s).\n")

  # a single (fully-directed) extension is returned as a bn object, several as a
  # list, matching the behaviour of cextend.all.backend().
  if (length(extensions) == 1)
    return(extensions[[1]])
  else
    return(extensions)

}#CEXTEND.ALL.PDAG.BACKEND

# maximally orient a PDAG into its MPDAG by applying Meek's four rules until no
# more edges can be oriented, keeping the directed arcs that are already present
# (Wienobst et al. 2023, Algorithm 7, line 4; Meek, 1995). The adjacency matrix
# is coded as amat[i, j] == 1 for an arc i -> j, so an undirected edge i - j is
# amat[i, j] == amat[j, i] == 1 and a directed arc i -> j is amat[i, j] == 1 with
# amat[j, i] == 0. Orienting i -> j means clearing the reverse entry amat[j, i].
meek.rules.closure = function(amat, nodes, debug = FALSE) {

  n = length(nodes)

  directed = function(i, j) (amat[i, j] == 1L) && (amat[j, i] == 0L)
  undirected = function(i, j) (amat[i, j] == 1L) && (amat[j, i] == 1L)
  adjacent = function(i, j) (amat[i, j] == 1L) || (amat[j, i] == 1L)

  repeat {

    changed = FALSE

    # Meek rule 1: a -> b - c with a and c non-adjacent implies b -> c.
    for (b in seq_len(n))
      for (cc in seq_len(n)) {

        if (!undirected(b, cc))
          next

        for (a in seq_len(n))
          if (directed(a, b) && (a != cc) && !adjacent(a, cc)) {

            amat[cc, b] = 0L
            changed = TRUE
            break

          }#THEN

      }#FOR

    # Meek rule 2: a -> b -> c with a - c implies a -> c.
    for (a in seq_len(n))
      for (cc in seq_len(n)) {

        if (!undirected(a, cc))
          next

        for (b in seq_len(n))
          if (directed(a, b) && directed(b, cc)) {

            amat[cc, a] = 0L
            changed = TRUE
            break

          }#THEN

      }#FOR

    # Meek rule 3: a - b, a - c, a - d with b -> c and d -> c, and b, d
    # non-adjacent, imply a -> c.
    for (a in seq_len(n))
      for (cc in seq_len(n)) {

        if (!undirected(a, cc))
          next

        neigh = which(sapply(seq_len(n), function(b)
                  (b != a) && (b != cc) && undirected(a, b) && directed(b, cc)))

        oriented = FALSE
        for (b in neigh) {
          for (d in neigh)
            if ((d != b) && !adjacent(b, d)) {

              amat[cc, a] = 0L
              changed = oriented = TRUE
              break

            }#THEN
          if (oriented)
            break
        }#FOR

      }#FOR

    # Meek rule 4: a - b, a - c, a - d with d -> c and c -> b, and b, d
    # non-adjacent, imply a -> b.
    for (a in seq_len(n))
      for (b in seq_len(n)) {

        if (!undirected(a, b))
          next

        oriented = FALSE
        for (cc in seq_len(n)) {

          if (!(undirected(a, cc) && directed(cc, b)))
            next

          for (d in seq_len(n))
            if ((d != b) && (d != a) && undirected(a, d) &&
                directed(d, cc) && !adjacent(b, d)) {

              amat[b, a] = 0L
              changed = oriented = TRUE
              break

            }#THEN

          if (oriented)
            break

        }#FOR

      }#FOR

    if (!changed)
      break

  }#REPEAT

  return(amat)

}#MEEK.RULES.CLOSURE

# split a maximally-oriented PDAG (MPDAG) into its buckets: the connected
# components of the undirected subgraph, each induced on the MPDAG so that the
# within-bucket directed arcs are retained. The directed arcs are classified as
# within-bucket (both endpoints in the same bucket, to be oriented consistently
# by the AMO enumeration) or between-bucket (endpoints in different buckets, to
# be re-attached verbatim to every extension). (Wienobst et al. 2023,
# Definition 1; Algorithm 7, lines 6-7.)
pdag.buckets = function(amat, nodes, debug = FALSE) {

  n = length(nodes)

  # undirected adjacency: an edge is present iff both directions are set.
  und = (amat == 1L) & (t(amat) == 1L)

  # label the connected components of the undirected subgraph by breadth-first
  # search, so that every node (isolated ones included) gets a bucket id.
  comp = integer(n)
  k = 0L
  for (s in seq_len(n)) {

    if (comp[s] != 0L)
      next

    k = k + 1L
    comp[s] = k
    queue = s

    while (length(queue) > 0) {

      v = queue[1]
      queue = queue[-1]
      reached = which(und[v, ] & (comp == 0L))
      comp[reached] = k
      queue = c(queue, reached)

    }#WHILE

  }#FOR

  names(comp) = nodes

  # the vertex set of each bucket.
  bucket.nodes = lapply(seq_len(k), function(i) nodes[comp == i])

  # classify the directed arcs by whether their endpoints share a bucket.
  dir.idx = which((amat == 1L) & (t(amat) == 0L), arr.ind = TRUE)
  same = comp[dir.idx[, 1]] == comp[dir.idx[, 2]]
  to.arcs = function(idx)
    matrix(c(nodes[idx[, 1]], nodes[idx[, 2]]), ncol = 2,
      dimnames = list(NULL, c("from", "to")))
  within.arcs = to.arcs(dir.idx[same, , drop = FALSE])
  between.arcs = to.arcs(dir.idx[!same, , drop = FALSE])

  # the undirected edges, each listed once (from < to).
  und.idx = which(und & upper.tri(und), arr.ind = TRUE)
  undirected.arcs = to.arcs(und.idx)

  if (debug) {

    cat("* found", k, "bucket(s):\n")
    for (i in seq_len(k))
      cat("  >", bucket.nodes[[i]], "\n")
    cat("* within-bucket directed arcs:", nrow(within.arcs), "\n")
    cat("* between-bucket directed arcs:", nrow(between.arcs), "\n")

  }#THEN

  return(list(components = comp, bucket.nodes = bucket.nodes,
    within.arcs = within.arcs, between.arcs = between.arcs,
    undirected.arcs = undirected.arcs))

}#PDAG.BUCKETS

# enumerate the AMOs of the buckets of an MPDAG with the S+ restriction: a
# modified Maximum Cardinality Search that, at each step, only picks a vertex from
# the highest-label set that has no unvisited within-bucket parent (the set S+).
# This makes every ordering conform to the within-bucket directed arcs, so that
# orienting the (undirected and directed) within-bucket edges by the ordering
# reproduces the directed arcs and lists each AMO exactly once. It is deliberately
# separate from uccg.all.extensions(), the plain (CPDAG) MCS enumerator.
# (Wienobst et al. 2023, Algorithm 6; Lemmas 3-4.) Returns a list of arc sets,
# each the within-bucket edges oriented by one valid ordering.
pdag.buckets.amos = function(buckets, nodes, debug = FALSE) {

  empty.arcs = matrix(character(0), ncol = 2,
                 dimnames = list(NULL, c("from", "to")))

  # the skeleton of the buckets: the undirected edges plus the within-bucket
  # directed arcs (with their direction dropped for the traversal, restored via
  # the ordering). every such arc lies inside a single bucket by construction.
  edges = rbind(buckets$undirected.arcs, buckets$within.arcs)

  # only the vertices carrying a within-bucket edge need to be oriented; the rest
  # (isolated vertices and singleton buckets) contribute a single empty AMO.
  active = intersect(nodes, unique(as.vector(edges)))
  if (length(active) == 0)
    return(list(empty.arcs))

  # neighbours in the skeleton and within-bucket parents of each active vertex.
  nbr = sapply(active, function(v) character(0), simplify = FALSE)
  for (r in seq_len(nrow(edges))) {

    nbr[[edges[r, 1]]] = c(nbr[[edges[r, 1]]], edges[r, 2])
    nbr[[edges[r, 2]]] = c(nbr[[edges[r, 2]]], edges[r, 1])

  }#FOR
  nbr = lapply(nbr, unique)

  parents = sapply(active, function(v) character(0), simplify = FALSE)
  for (r in seq_len(nrow(buckets$within.arcs)))
    parents[[buckets$within.arcs[r, 2]]] =
      c(parents[[buckets$within.arcs[r, 2]]], buckets$within.arcs[r, 1])

  # orient every skeleton edge by the position of its endpoints in the ordering.
  orient.by = function(ordering) {

    pos = structure(seq_along(ordering), names = ordering)
    oriented = edges
    swap = pos[edges[, 1]] > pos[edges[, 2]]
    oriented[swap, ] = edges[swap, c(2, 1)]
    dimnames(oriented) = list(NULL, c("from", "to"))

    return(oriented)

  }#ORIENT.BY

  # vertices reachable from 'start' using only skeleton edges among 'allowed'.
  reachable = function(start, allowed) {

    seen = start
    queue = start

    while (length(queue) > 0) {

      v = queue[1]
      queue = queue[-1]
      fresh = setdiff(intersect(nbr[[v]], allowed), seen)
      seen = c(seen, fresh)
      queue = c(queue, fresh)

    }#WHILE

    return(setdiff(seen, start))

  }#REACHABLE

  # the recursion state and the results, held by reference in an environment.
  n = length(active)
  state = new.env()
  state$A = vector(n + 1, mode = "list")
  state$A[[1]] = active
  state$tau = character(0)
  state$results = vector(0, mode = "list")

  # locate the label set that currently contains vertex w.
  index.of = function(w)
    which(vapply(state$A, function(s) w %in% s, logical(1)))

  recurse = function() {

    # a complete ordering: orient the within-bucket edges and store the AMO.
    if (length(state$tau) == n) {

      state$results = c(state$results, list(orient.by(state$tau)))
      return(invisible(NULL))

    }#THEN

    # the highest-label set, restricted to the vertices with no unvisited parent.
    i = max(which(lengths(state$A) > 0))
    Splus = state$A[[i]][vapply(state$A[[i]],
              function(v) all(parents[[v]] %in% state$tau), logical(1))]

    if (length(Splus) == 0)
      stop("no S+ vertex is available: 'x' is not an extendable MPDAG.")

    v = Splus[1]
    x = v
    R = character(0)

    repeat {

      # add x to the ordering and raise the label of its unvisited neighbours.
      state$A[[i]] = setdiff(state$A[[i]], x)
      state$tau = c(state$tau, x)

      if (debug)
        cat("  @ considering ordering:", state$tau, "\n")

      for (w in setdiff(nbr[[x]], state$tau)) {

        j = index.of(w)
        state$A[[j]] = setdiff(state$A[[j]], w)
        state$A[[j + 1]] = c(state$A[[j + 1]], w)

      }#FOR

      recurse()

      # undo: lower the labels back and remove x from the ordering.
      for (w in setdiff(nbr[[x]], state$tau)) {

        j = index.of(w)
        state$A[[j]] = setdiff(state$A[[j]], w)
        state$A[[j - 1]] = c(state$A[[j - 1]], w)

      }#FOR

      state$A[[i]] = c(state$A[[i]], x)
      state$tau = setdiff(state$tau, x)

      # explore the other first-vertex choices reachable from v within S+ only
      # (choices outside S+ or unreachable from v would duplicate AMOs).
      if (x == v)
        R = reachable(v, Splus)

      if (length(R) == 0)
        break

      x = R[length(R)]
      R = R[-length(R)]

    }#REPEAT

  }#RECURSE

  recurse()

  return(state$results)

}#PDAG.BUCKETS.AMOS

# enumerate all possible extensions of an undirected connected chordal graph.
uccg.all.extensions = function (ug, state, generate = TRUE, debug = FALSE) {

  n = nnodes(ug)

  # initialise the global state at the top of the recursion.
  if (missing(state)) {

    state = new.env()
    state$A = vector(n, mode = "list")
    state$A[[1]] = nodes(ug)

    state$tau = character(0)

    state$count = 0
    state$AMOs = vector(0, mode = "list")

  }#THEN

  # tau contains a complete ordering, produce the directed extension.
  if (length(state$tau) == n) {

    state$count = state$count + 1

    if (debug)
      cat("  @ found consistent extension #", state$count, "!\n")

    if (generate) {

      state$AMOs =
        c(state$AMOs, list(pdag2dag.backend(ug$arcs, ordering = state$tau)))

    }#THEN

    return(NULL)

  }#THEN

  # here starts the actual recursive exploration.
  i = max(which(sapply(state$A, length) > 0))
  v = state$A[[i]][1]
  x = v

  repeat {

    # add node x to the node ordering, and remove it from the candidates.
    state$A[[i]] = setdiff(state$A[[i]], x)
    state$tau = c(state$tau, x)

    if (debug)
      cat("* considering ordering:", state$tau, "\n")

    for (w in setdiff(ug$nodes[[x]]$nbr, state$tau)) {

      j = which(sapply(state$A, `%in%`, x = w))
      state$A[[j]] = setdiff(state$A[[j]], w)
      state$A[[j + 1]] = c(state$A[[j + 1]], w)

    }#FOR

    # recurse.
    uccg.all.extensions(ug, state, generate = generate, debug = debug)

    # remove node x from the node ordering.
    for (w in setdiff(ug$nodes[[x]]$nbr, state$tau)) {

      j = which(sapply(state$A, `%in%`, x = w))
      state$A[[j]] = setdiff(state$A[[j]], w)
      state$A[[j - 1]] = c(state$A[[j - 1]], w)

    }#FOR

    state$A[[i]] = c(state$A[[i]], x)
    state$tau = setdiff(state$tau, x)

    # update the node search space.
    if (x == v) {

      subg = subgraph.backend(ug, state$A[[i]])
      R = sapply(setdiff(names(subg$nodes), v), path.exists, x = subg, from = v)

      if (length(R) == 0)
        R = character(0)
      else
        R = names(which(R))

    }#THEN

    if (length(R) == 0)
      break
    else {

      x = R[length(R)]
      R = R[-length(R)]

    }#ELSE

  }#REPEAT

  if (generate)
    return(state$AMOs)
  else
    return(state$count)

}#UCCG.ALL.EXTENSIONS
