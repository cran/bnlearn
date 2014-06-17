
# dump a bn.fit object into a BIF/DSC/NET file.
write.foreign.backend = function(fd, fitted, format = "bif") {

  # print the preamble.
  if (format == "dsc")
    cat("belief network \"unknown\"\n", file = fd)
  if (format == "net")
    cat("net \n{ \n}\n", file = fd)
  else if (format == "bif")
    cat("network unknown {\n}\n", file = fd)

  # get the levels and the parents of each node.
  levels = sapply(names(fitted), function(x) dimnames(fitted[[x]]$prob)[[1]],
             simplify = FALSE)
  levels = sanitize.levels(levels)

  parents = lapply(fitted, "[[", "parents")

  # print the variable decalarations, describing the number and the labels of
  # the levels of eah node.
  for (node in names(fitted)) {

    if (format == "dsc")
      dsc.write.node(node = node, levels = levels[[node]], fd = fd)
    else if (format == "net")
      net.write.node(node = node, levels = levels[[node]], fd = fd)
    else if (format == "bif")
      bif.write.node(node = node, levels = levels[[node]], fd = fd)

  }#FOR

  # print the (conditional) probbability table of each node.
  for (node in names(fitted)) {

    cpt = fitted[[node]]$prob
    nlevels = length(levels[[node]])

    if (format == "dsc")
      dsc.write.probabilities(node = node, nlevels = nlevels, cpt = cpt,
        parents = parents[[node]], fd = fd)
    else if (format == "net")
      net.write.probabilities(node = node, nlevels = nlevels, cpt = cpt,
        parents = parents[[node]], fd = fd)
    else if (format == "bif")
      bif.write.probabilities(node = node, nlevels = nlevels, cpt = cpt,
        parents = parents[[node]], fd = fd)

  }#FOR

}#WRITE.BIF.BACKEND

sanitize.levels = function(levels) {

  clean = lapply(levels, gsub, pattern = "[,{}()#%]", replacement = "_")

  # preserve names, which is useful when the input comes from dimnames().
  names(clean) = names(levels)

  if (any(sapply(clean, anyDuplicated) > 0))
    stop("duplicated levels after sanitization.")

  return(clean)

}#SANITIZE.LEVELS

bif.write.node = function(node, levels, fd) {

  cat("variable", node, "{\n  type discrete [", length(levels),
    "] {", paste(levels, collapse = ", "),
    "};\n}\n", file = fd)

}#BIF.WRITE.NODE

bif.write.probabilities = function(node, nlevels, cpt, parents, fd) {

  if (length(parents) == 0) {

    probs = paste(format(cpt, nsmall = 1), collapse = ", ")

    cat("probability (", node, ") {\n", file = fd)
    cat("  table", paste(probs, ";\n", sep = ""), file = fd)
    cat("}\n", file = fd)


  }#THEN
  else {

    configs = sanitize.levels(dimnames(cpt)[-1])
    configs = expand.grid(configs, stringsAsFactors = FALSE)
    configs = apply(configs, 1, function(x) paste("(",
                paste(x, collapse = ", "), ")", sep = "") )

    cat("probability (", node, "|", paste(parents, collapse = ", "),
      ") {\n", file = fd)

    for (i in seq_along(configs)) {

      probs = cpt[nlevels * (i - 1) + seq(nlevels)]

      # paper over missing values with a uniform distribution
      if (all(is.na(probs)))
        probs = rep(1, length(probs))/length(probs)

      probs = paste(format(probs, nsmall = 1), collapse = ", ")
      cat(" ", configs[i], paste(probs, ";\n", sep = ""), file = fd)

    }#FOR

    cat("}\n", file = fd)

  }#THEN

}#BIF.WRITE.PROBABILITIES

net.write.node = function(node, levels, fd) {

  cat("node", node, "\n{\n  states = (", paste("\"", paste(levels,
        collapse = "\" \""), "\"", sep = ""), ");\n}\n", file = fd)

}#NET.WRITE.NODE

net.write.probabilities = function(node, nlevels, cpt, parents, fd) {

  if (length(parents) == 0) {

    probs = paste(format(cpt, nsmall = 1))

    cat("potential (", node, ") \n{\n", file = fd)
    cat("  data = (", probs, ");\n}\n", file = fd)

  }#THEN
  else {

    cat("potential (", node, "|", paste(parents, collapse = " "),
      ") \n{\n  data = ", file = fd)

    # in NET files the "most significant" variable for the ordering is the
    # last one, not the first one.
    configs = sanitize.levels(dimnames(cpt))
    configs = expand.grid(configs[1+rev(1:length(parents))])
    configs = configs[rev(1:ncol(configs))]
    indexes = configurations(configs, factor = FALSE)

    # paste CPT columns together and save them in the right order.
    probs = character(length(indexes))

    for (i in seq_along(indexes)) {

      temp.probs = cpt[nlevels * (indexes[i] - 1) + seq(nlevels)]

      # paper over missing values with a uniform distribution
      if (all(is.na(temp.probs)))
        temp.probs = rep(1, length(temp.probs))/length(temp.probs)

      probs[i] = paste(format(temp.probs, nsmall = 1), collapse = " ")
      probs[i] = paste("(", probs[i], ")", sep = "")

    }#FOR

    # add braces in the proper places.
    parents.nlevels = dim(cpt)[-1]

    for (i in rev(parents.nlevels)) {

      blocks = (i - 1 + seq_along(probs)) %/% i
      blocks = lapply(seq(max(blocks)), function(x) probs[blocks == x])
      probs = sapply(blocks, function(x) paste("(", paste(x, collapse = ""), ")",
                sep = "", collapse = ""))

    }#FOR

    cat(probs, ";\n}\n", file = fd)

  }#ELSE

}#NET.WRITE.PROBABILITIES

dsc.write.node = function(node, levels, fd) {

  cat("node", node, "{\n  type : discrete [", length(levels),
    "] = {", paste("\"", paste(levels, collapse = "\", \""), "\"", sep = ""),
    "};\n}\n", file = fd)

}#DSC.WRITE.NODE

dsc.write.probabilities = function(node, nlevels, cpt, parents, fd) {

  if (length(parents) == 0) {

    probs = paste(format(cpt, nsmall = 1), collapse = ", ")

    cat("probability (", node, ") {\n", file = fd)
    cat("  ", paste(probs, ";\n", sep = ""), file = fd)
    cat("}\n", file = fd)

  }#THEN
  else {

    configs = expand.grid(lapply(dim(cpt)[-1], function(x) seq(x) - 1), stringsAsFactors = FALSE)
    configs = apply(configs, 1, function(x) paste("(", paste(x, collapse = ", "), ")", sep = "") )

    cat("probability (", node, "|", paste(parents, collapse = ", "),
      ") {\n", file = fd)

    for (i in seq_along(configs)) {

      probs = cpt[nlevels * (i - 1) + seq(nlevels)]

      # paper over missing values with a uniform distribution
      if (all(is.na(probs)))
        probs = rep(1, length(probs))/length(probs)

      probs = paste(format(probs, nsmall = 1), collapse = ", ")
      cat(" ", configs[i], ":", paste(probs, ";\n", sep = ""), file = fd)

    }#FOR

    cat("}\n", file = fd)

  }#ELSE

}#DSC.WRITE.PROBABILITIES

# export graph structures as DOT files.
write.dot.backend = function(fd, graph) {

  if (is(graph, "bn"))
    arcs = graph$arcs
  else
    arcs = fit2arcs(graph)

  nodes = .nodes(graph)

  if (nrow(arcs) == 0) {

    # this is an empty graph, we assume it's (partially) directed.
    cat("digraph {\n", file = fd)
    cat(paste('  "', nodes(graph), '" ;'), sep = "\n")
    cat("}\n", file = fd)

  }#THEN
  else {

    wund = which.undirected(arcs, nodes)

    type = ifelse(length(which(wund)) == nrow(arcs), "graph", "digraph")

    if (type == "graph") {

      dedup = match(arcs[, "from"], nodes) < match(arcs[, "to"], nodes)

      und = apply(arcs[wund & dedup, , drop = FALSE], 1,
              function(x) paste('"', x, '"', sep = "", collapse = " -- "))

      # this is an undirected graph.
      cat("graph {\n", file = fd)
      cat(paste('  "', nodes(graph), '" ;', sep = ""), sep = "\n", file = fd)
      cat(paste("  edge [dir=none]", und, ";\n"), sep = "", file = fd)
      cat("}\n", file = fd)

    }#THEN
    else if (type == "digraph") {

      # this is a (partially) directed graph.
      cat("digraph {\n", file = fd)

      dir = apply(arcs[!wund, , drop = FALSE], 1,
              function(x) paste('"', x, '"', sep = "", collapse = " -> "))
      cat(paste('  "', nodes(graph), '" ;', sep = ""), sep = "\n", file = fd)

      if (any(wund)) {

        dedup = match(arcs[, "from"], nodes) < match(arcs[, "to"], nodes)
        und = apply(arcs[wund & dedup, , drop = FALSE], 1,
                function(x) paste('"', x, '"', sep = "", collapse = " -> "))
        cat(paste("  edge [dir=none]", und, ";\n"), sep = "", file = fd)

      }#THEN

      cat(paste("  edge [dir=forward]", dir, ";\n"), sep = "", file = fd)
      cat("}\n", file = fd)

    }#ELSE

  }#ELSE

}#WRITE.DOT.BACKEND

