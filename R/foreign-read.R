
# read a BIF file into a bn.fit object.
read.foreign.backend = function(lines, format = "bif", filename, debug = FALSE) {

  # remove comments.
  if (format %in% c("bif", "dsc"))
    lines = sub("//.*", "", lines)
  else if (format == "net")
    lines = sub("%.*", "", lines)
  # remove indentation and whitespace damage.
  lines = gsub("^\\s+|\\s+$", "", lines, perl = TRUE)
  # remove empty lines.
  lines = lines[grep("^\\s*$", lines, perl = TRUE, invert = TRUE)]

  if (format == "bif") {

    # remove useless properties from declarations and CPTs.
    lines = bif.preparse(lines)

  }#THEN
  else if (format == "dsc") {

    # quick fix for CPT declaration splitting (thanks Genie).
    lines = dsc.preparse(lines)

  }#THEN

  # check the file banner.
  if (format == "bif")
    bif.check.banner(banner = lines[1], filename = filename)
  else if (format == "dsc")
    dsc.check.banner(banner = lines[1], filename = filename)
  else if (format == "net")
    net.check.banner(banner = lines[1], filename = filename)

  # get the node labels and the size of node set.
  if (format == "bif")
    nodes = bif.get.nodes(lines)
  else if (format == "dsc")
    nodes = dsc.get.nodes(lines)
  else if (format == "net")
    nodes = net.get.nodes(lines)
  nnodes = length(nodes)

  # check whether all variables are discrete.
  if (format == "bif")
    bif.check.discrete(lines, nnodes)
  else if (format == "dsc")
    dsc.check.discrete(lines, nnodes)
  else if (format == "net")
    net.check.discrete(lines, nnodes)

  # check whether the labels on the probability tables match with the node labels.
  if (format == "bif")
    cpts = bif.get.cpt.names(lines)
  else if (format == "dsc")
    cpts = dsc.get.cpt.names(lines)
  else if (format == "net")
    cpts = net.get.cpt.names(lines)

  if (!setequal(cpts, nodes)) {

    missing.cpts = setdiff(nodes, cpts)
    missing.nodes = setdiff(cpts, nodes)
    bogus = unique(c(missing.cpts, missing.nodes))

    for (m in missing.cpts)
      warning("the CPT corresponding to node ", m, " is missing, dropping.")

    for (m in missing.nodes)
      warning("the node description of node ", m, " is missing, dropping.")

    # recompute fundamental quantities.
    nodes = nodes[!(nodes %in% bogus)]
    nnodes = length(nodes)

  }#THEN
  else {

    bogus = character(0)

  }#ELSE

  # find out where each node description begins.
  if (format == "bif")
    description.start = bif.get.node.descriptions(lines, nodes)
  else if (format == "dsc")
    description.start = dsc.get.node.descriptions(lines, nodes)
  else if (format == "net")
    description.start = net.get.node.descriptions(lines, nodes)

  # find out where each conditional probability table begins.
  if (format == "bif")
    cpt.start = bif.get.cpt.descriptions(lines, cpts)
  else if (format == "dsc")
    cpt.start = dsc.get.cpt.descriptions(lines, cpts)
  else if (format == "net")
    cpt.start = net.get.cpt.descriptions(lines, cpts)

  # get the levels associated with each node.
  if (format == "bif") {

    nodes.levels = sapply(nodes, bif.get.levels, start = description.start,
                     lines = lines, simplify = FALSE)

  }#THEN
  else if (format == "dsc") {

    nodes.levels = sapply(nodes, dsc.get.levels, start = description.start,
                     lines = lines, simplify = FALSE)

  }#THEN
  else if (format == "net") {

    nodes.levels = sapply(nodes, net.get.levels, start = description.start,
                     lines = lines, simplify = FALSE)

  }#THEN

  # all nodes should have at least two levels, drop dummy nodes with a warning.
  dummies = names(which(sapply(nodes.levels, length) < 2))

  if (length(dummies) > 0) {

    for (d in dummies)
      warning("node ", d, " have only one level, dropping.")

    # recompute some fundamental quantities.
    nodes = nodes[!(nodes %in% dummies)]
    nnodes = length(nodes)
    description.start = description.start[nodes]
    cpt.start = cpt.start[nodes]

  }#THEN

  # get the parents of each node.
  if (format == "bif") {

    parents = bif.get.parents(lines, start = cpt.start, dummies = dummies,
                bogus = bogus)

  }#THEN
  else if (format == "dsc") {

    parents = dsc.get.parents(lines, start = cpt.start, dummies = dummies,
                bogus = bogus)
  }#THEN
  else if (format == "net") {

    parents = net.get.parents(lines, start = cpt.start, dummies = dummies,
                bogus = bogus)
  }#THEN

  # separate root and non-root nodes.
  nonroot.nodes = names(parents)
  root.nodes = nodes[!(nodes %in% nonroot.nodes)]

  # create the empty bn.fit object.
  fitted = structure(vector(nnodes, mode = "list"), names = nodes)

  # fill in the metadata for the root nodes.
  for (node in nodes) {

    # get the parent set.
    if (node %in% root.nodes)
      parent.set = character(0)
    else
      parent.set = parents[[node]]

    if (debug) {

      if (node %in% root.nodes)
        cat("* found root node", node, ".\n")
      else
        cat("* found node", node, "with parents", parent.set, ".\n")
      cat("  > node", node, "has levels", nodes.levels[[node]], "\n")

    }#THEN

    # build the (conditional) probability table.
    if (format == "bif") {

      node.cpt = bif.get.probabilities(node, start = cpt.start, lines = lines,
                   nodes.levels = nodes.levels, parents = parents,
                   root = (node %in% root.nodes))

    }#THEN
    else if (format == "dsc") {

      node.cpt = dsc.get.probabilities(node, start = cpt.start, lines = lines,
                   nodes.levels = nodes.levels, parents = parents,
                   root = (node %in% root.nodes))

    }#THEN
    else if (format == "net") {

      node.cpt = net.get.probabilities(node, start = cpt.start, lines = lines,
                   nodes.levels = nodes.levels, parents = parents,
                   root = (node %in% root.nodes))

    }#THEN

    fitted[[node]] = structure(list(node = node, parents = parent.set,
                       children = foreign.get.children(node, parents),
                       prob = node.cpt), class = "bn.fit.dnode")

    if (debug) {

      cat("  > conditional probability table:\n")
      print(node.cpt)

    }#THEN

  }#FOR

  # set the class of the return value; doing it before the object has been
  # completely initialized clashed with the "[[<-.bn.fit" replacement methods.
  class(fitted) = "bn.fit"

  return(fitted)

}#READ.FOREIGN.BACKEND

# remove properties to avoid confusion.
bif.preparse = function(lines) {

  # remove properties.
  lines = gsub("^\\s*property[^;]*;", "", lines)
  # remove empty lines.
  lines = lines[grep("^\\s*$", lines, perl = TRUE, invert = TRUE)]

  return(lines)

}#BIF.PREPARSE

# rebuild declarations split over multiple lines.
dsc.preparse = function(lines) {

  p = grep("probability[^)]*$", lines)
  for (l in p) {

    end = match.brace(lines, start = l, open = "(", close = ")")
    lines[l] = paste(lines[l:end], collapse = " ")
    lines[l] = gsub(" , ", ", ", lines[l])
    lines[(l+1):end] = ""

  }
  # remove empty lines.
  lines = lines[grep("^\\s*$", lines, perl = TRUE, invert = TRUE)]

  return(lines)

}#DSC.PREPARSE

# check the banner of a BIF file.
bif.check.banner = function(banner, filename) {

  if (length(grep("^network", banner)) != 1)
    stop("the file '", filename, "' does not conform to the BIF format.")

}#BIF.CHECK.BANNER

# check the banner of a DSC file.
dsc.check.banner = function(banner, filename) {

  if (length(grep("^belief network", banner)) != 1)
    stop("the file '", filename, "' does not conform to the DSC format.")

}#DSC.CHECK.BANNER

# check the banner of a NET file.
net.check.banner = function(banner, filename) {

  if (length(grep("^net", banner)) != 1)
    stop("the file '", filename, "' does not conform to the NET format.")

}#NET.CHECK.BANNER

# get the node labels from a BIF file.
bif.get.nodes = function(lines) {

  sub("variable\\s+(.+)\\s+\\{", "\\1",
    grep("^variable .+", lines, value = TRUE), perl = TRUE)

}#BIF.GET.NODES

# get the node labels from a DSC file.
dsc.get.nodes = function(lines) {

  sub("node\\s+(\\w+)\\s*\\{*", "\\1",
    grep("^node .+", lines, value = TRUE), perl = TRUE)

}#DSC.GET.NODES

# get the node labels from a NET file.
net.get.nodes = function(lines) {

  sub("node\\s+(\\w+)\\s*\\{*", "\\1",
    grep("^node .+", lines, value = TRUE), perl = TRUE)

}#NET.GET.NODES

# count the number of discrete nodes in a BIF file.
bif.check.discrete = function(lines, nnodes) {

  ndisc = length(grep("type\\s+discrete", lines))

  if (ndisc != nnodes)
    stop("only BIF files describing discrete networks are supported.")

}#BIF.CHECK.DISCRETE

# count the number of discrete nodes in a DSC file.
dsc.check.discrete = function(lines, nnodes) {

  ndisc = length(grep("type\\s*[:]{0,1}\\s*discrete", lines))

  if (ndisc != nnodes)
    stop("only DSC files describing discrete networks are supported.")

}#DSC.CHECK.DISCRETE

# count the number of discrete nodes in a DSC file.
net.check.discrete = function(lines, nnodes) {

  if (length(grep("^continuous\\s+node\\s+\\w+", lines)) > 0)
    stop("only NET files describing discrete networks are supported.")
  if (length(grep("^decision\\s+\\w+", lines)) > 0)
    stop("NET files describing decision networks are not supported.")

}#NET.CHECK.DISCRETE

# get the labels of the CPTs from a BIF file, to compare them to node labels.
bif.get.cpt.names = function(lines) {

  sub("probability\\s*\\(\\s*(\\w+).*", "\\1",
      grep("^probability", lines, value = TRUE))

}#BIF.GET.CPT.NAMES

# get the labels of the CPTs from a DSC file, to compare them to node labels.
dsc.get.cpt.names = bif.get.cpt.names

# get the labels of the CPTs from a NET file, to compare them to node labels.
net.get.cpt.names = function(lines) {

  sub("potential\\s*\\(\\s*(\\w+).*", "\\1",
      grep("^potential", lines, value = TRUE))

}#NET.GET.CPT.NAMES

# get the starting lines of node descriptions in a BIF file.
bif.get.node.descriptions = function(lines, nodes) {

  desc = grep(paste("variable\\s+(", paste(nodes, collapse = "|"), ")", sep = ""), lines)
  names(desc) = nodes

  return(desc)

}#BIF.GET.NODE.DESCRIPTIONS

# get the starting lines of node descriptions in a DSC file.
dsc.get.node.descriptions = function(lines, nodes) {

  desc = grep(paste("node\\s+(", paste(nodes, collapse = "|"), ")", sep = ""), lines)
  names(desc) = nodes

  return(desc)

}#DSC.GET.NODE.DESCRIPTIONS

# get the starting lines of node descriptions in a NET file.
net.get.node.descriptions = function(lines, nodes) {

  desc = grep(paste("node\\s+(", paste(nodes, collapse = "|"), ")", sep = ""), lines)
  names(desc) = nodes

  return(desc)

}#NET.GET.NODE.DESCRIPTIONS

# get the starting lines of CPTs in a BIF file.
bif.get.cpt.descriptions = function(lines, cpts) {

  desc = grep(paste("probability\\s*\\(\\s*(", paste(cpts, collapse = "|"), ")", sep = ""), lines)
  names(desc) = cpts

  return(desc)

}#BIF.GET.CPT.DESCRIPTIONS

# get the starting lines of CPTs in a DSC file.
dsc.get.cpt.descriptions = bif.get.cpt.descriptions

# get the starting lines of CPTs in a NET file.
net.get.cpt.descriptions = function(lines, cpts) {

  desc = grep(paste("potential\\s*\\(\\s*(", paste(cpts, collapse = "|"), ")", sep = ""), lines)
  names(desc) = cpts

  return(desc)

}#NET.GET.CPT.DESCRIPTIONS

# get the parents of each node in a BIF file.
bif.get.parents = function(lines, start, dummies, bogus) {

  # get the dependencies.
  parents = (lines[start])[grep("\\|", lines[start])]
  # extract the node labels and the labels of the respective parents.
  parents = sub("probability\\s+\\(\\s+(.+)\\s+\\|(.+)\\)\\s+\\{", "\\1 \\2",
               parents)
  # split the labels of the parents.
  nonroot.nodes = sapply(strsplit(parents, " "), "[", 1)
  parents = sapply(strsplit(parents, " "), "[", -1, simplify = FALSE)
  names(parents) = nonroot.nodes
  # remove commas and empty values.
  parents = lapply(parents, function(x) {

    p = sub(",", "", x)
    p = p[grep("^\\s*$", p, perl = TRUE, invert = TRUE)]

    # check whether there are dropped nodes among the parents.
    if (any(p %in% dummies))
      stop("dropped dummy node is the parent of another node.")
    if (any(p %in% bogus))
      stop("dropped mismatched node is the parent of another node.")

    return(p)

  } )

  return(parents)

}#BIF.GET.PARENTS

# get the parents of each node in a DSC file.
dsc.get.parents = function(lines, start, dummies, bogus) {

  # get the dependencies.
  parents = (lines[start])[grep("\\|", lines[start])]
  # extract the node labels and the labels of the respective parents.
  parents = sub("probability\\s*\\(\\s*(.+)\\s*\\|\\s*(.+)\\s*\\).*", "\\1 \\2",
               parents)
  # split the labels of the parents.
  nonroot.nodes = sapply(strsplit(parents, " "), "[", 1)
  parents = sapply(strsplit(parents, " "), "[", -1, simplify = FALSE)
  names(parents) = nonroot.nodes
  # remove commas and empty values.
  parents = lapply(parents, function(x) {

    p = sub(",", "", x)
    p = p[grep("^\\s*$", p, perl = TRUE, invert = TRUE)]

    # check whether there are dropped nodes among the parents.
    if (any(p %in% dummies))
      stop("dropped dummy node is the parent of another node.")
    if (any(p %in% bogus))
      stop("dropped mismatched node is the parent of another node.")

    return(p)

  } )

  return(parents)

}#DSC.GET.PARENTS

# get the parents of each node in a NET file.
net.get.parents = function(lines, start, dummies, bogus) {

  # get the dependencies (robust against Genie output).
  parents = (lines[start])[grep("\\|\\s*\\w+", lines[start])]
  # extract the node labels and the labels of the respective parents.
  parents = sub("potential\\s+\\(\\s*(.+)\\s*\\|\\s*(.+)\\s*\\).*", "\\1 \\2",
               parents)
  # split the labels of the parents.
  nonroot.nodes = sapply(strsplit(parents, " "), "[", 1)
  parents = sapply(strsplit(parents, " "), "[", -1, simplify = FALSE)
  names(parents) = nonroot.nodes
  # remove commas and empty values.
  parents = lapply(parents, function(x) {

    p = sub(",", "", x)
    p = p[grep("^\\s*$", p, perl = TRUE, invert = TRUE)]

    # check whether there are dropped nodes among the parents.
    if (any(p %in% dummies))
      stop("dropped dummy node is the parent of another node.")
    if (any(p %in% bogus))
      stop("dropped mismatched node is the parent of another node.")

    return(p)

  } )

  return(parents)

}#NET.GET.PARENTS

foreign.get.children = function(node, parents) {

  names(which(sapply(parents, function(x, node) { node %in% x  }, node = node)))

}#BIF.GET.CHILDREN

# get the levels of each node in a BIF file.
bif.get.levels = function(node, start, lines) {

  end.line = 0

  # get the line in which the description is starting.
  start.line = start[node]
  # loop until the line in which the description is ending.
  end.line = match.brace(lines, start.line)

  # get all the node's description on one line for easy handling.
  desc = paste(lines[start.line:end.line], collapse = "")
  # deparse the node's level.
  levels = sub(".+type\\s+discrete\\s*\\[\\s*\\d+\\s*\\]\\s+[=]*\\s*\\{\\s+(.+)\\s*\\}.+", "\\1", desc)
  levels = strsplit(levels, ",")[[1]]
  levels = sub("^\\s*(.+?)\\s*$", "\\1", levels)

  if (any(duplicated(levels)))
    stop("duplicated levels '", paste(levels[duplicated(levels)], collapse = "' '"), "' for node ", node, ".\n")

  return(levels)

}#BIF.GET.LEVELS

# get the levels of each node in a DSC file.
dsc.get.levels = function(node, start, lines) {

  end.line = 0

  # get the line in which the description is starting.
  start.line = start[node]
  # loop until the line in which the description is ending.
  end.line = match.brace(lines, start.line)

  # get all the node's description on one line for easy handling.
  desc = paste(lines[start.line:end.line], collapse = "")
  # deparse the node's level.
  levels = sub(".+type\\s*[:]{0,1}\\s*discrete\\s*\\[\\s*\\d+\\s*\\]\\s+[=]{0,1}\\s*\\{\\s*(.+)\\s*\\}.+", "\\1", desc)
  levels = strsplit(levels, "\"\\s*,")[[1]]
  levels = sub("^\\s*\"*(.+?)\"*\\s*$", "\\1", levels)

  if (any(duplicated(levels)))
    stop("duplicated levels '", paste(levels[duplicated(levels)], collapse = "' '"), "' for node ", node, ".\n")

  return(levels)

}#DSC.GET.LEVELS

# get the levels of each node in a NET file.
net.get.levels = function(node, start, lines) {

  end.line = 0

  # get the line in which the description is starting.
  start.line = start[node]
  # loop until the line in which the description is ending.
  end.line = match.brace(lines, start.line)

  # get all the node's description on one line for easy handling.
  desc = paste(lines[start.line:end.line], collapse = "")

  # deparse the node's level.
  levels = sub(".+states\\s*[=]{0,1}\\s*\\(\\s*(.+?)\\s*\\).+", "\\1", desc)
  levels = gsub("\"", "", strsplit(levels, "\"\\s*\"")[[1]])

  if (any(duplicated(levels)))
    stop("duplicated levels '", paste(levels[duplicated(levels)], collapse = "' '"), "' for node ", node, ".\n")

  return(levels)

}#NET.GET.LEVELS

# build the (conditional) probability table for a node in a BIF file.
bif.get.probabilities = function(node, start, lines, nodes.levels, parents, root) {

  end.line = 0

  # get the line in which the description is starting.
  start.line = start[node]
  # loop until the line in which the description is ending.
  end.line = match.brace(lines, start.line)

  # get all the node's description on one line for easy handling.
  desc = paste(lines[start.line:end.line], collapse = "")

  # deparse the node's probability table.
  if (root) {

    probs = sub(".+table\\s+(.+)\\s*;\\s*\\}.*", "\\1", desc)
    probs = strsplit(probs, ",")[[1]]

    if (length(probs) != length(nodes.levels[[node]]))
      stop("the dimension of the CPT of node ", node,
        " does not match the number of its levels.")

    node.cpt = as.table(as.numeric(probs))
    dimnames(node.cpt) = list(nodes.levels[[node]])

  }#THEN
  else {

    row = strsplit(sub(".+\\{[^(]*(\\(.+?)\\s*[;]*\\s*\\}.*", "\\1", desc), ";")[[1]]

    # extract the conditional probability distributions.
    probs = strsplit(sub(".*\\)\\s*(.+)", "\\1", row), ",")
    probs = lapply(probs, as.numeric)

    # extract the configurations of the parents.
    cfg = strsplit(sub(".*\\((.+)\\).+", "\\1", row), ",")
    cfg = lapply(cfg, sub, pattern = "^\\s*(.+?)\\s*$", replacement = "\\1")
    dims = lapply(c(node, parents[[node]]), function(x) nodes.levels[[x]])

    # check whether the number of conditional probability distributions matches
    # the number of configurations.
    if (length(probs) != length(cfg))
      stop("the number of conditional distributions for node ", node,
        " do not math the number of configurations of its parents")
    # check whether the number of configurations matches the expected number of
    # configurations given the parents' number of levels.
    expected.cfgs = prod(sapply(nodes.levels[parents[[node]]], length))
    if (length(probs) != expected.cfgs)
      stop("there should be ", expected.cfgs, " conditional probability ",
        "distributions, but only ", length(cfg), " are present.")
    # check whether each conditional probability distribution is valid.
    not.prob = !sapply(probs, is.probability.vector)
    if (any(not.prob)) {

      not.prob = paste(which(not.prob), collapse = ", ")
      stop("conditional probability ditribution(s) ", not.prob, " of node ",
        node, " is not a vector of probabilities.")

    }#THEN
    # check whether each conditional probability distribution sums to one.
    not.prob = (abs(sapply(probs, sum) - 1) > 0.011)
    if (any(not.prob)) {

      # if more than 1% of probability mass, let's assume that the conditional
      # probability distribution is misspecied.
      not.prob = paste(which(not.prob), collapse = ", ")
      stop("conditional probability ditribution(s) ", not.prob, " of node ",
        node, " does not sum to one.")

    }#THEN
    else {

      # perform some more rounding to make the total probability mass closer
      # to one.
      probs = lapply(probs, prop.table)

    }#ELSE
    # check whether the conditional probability distributions have the right
    # number of elements
    if (any(lapply(probs, function(x) length(x)) != length(nodes.levels[[node]])))
      stop("one of the conditional probability ditributions of node ", node,
        " has the wrong number of elements.")

    node.cpt = table(seq(prod(sapply(dims, length))))
    dim(node.cpt) = sapply(dims, length)
    dimnames(node.cpt) = dims

    for (i in seq_along(probs)) {

      node.cpt = do.call("[<-", c(list(node.cpt, 1:length(nodes.levels[[node]])),
                   cfg[[i]], list(probs[[i]])))

    }#FOR

  }#ELSE

  return(node.cpt)

}#BIF.GET.PROBABILITIES

# build the (conditional) probability table for a node in a DSC file.
dsc.get.probabilities = function(node, start, lines, nodes.levels, parents, root) {

  end.line = 0

  nlevels = length(nodes.levels[[node]])

  # get the line in which the description is starting.
  start.line = start[node]
  # loop until the line in which the description is ending.
  end.line = match.brace(lines, start.line)

  # get all the node's description on one line for easy handling.
  desc = paste(lines[start.line:end.line], collapse = "")

  # deparse the node's probability table.
  if (root) {

    probs = sub(".+\\{\\s*(.+)\\s*;\\s*\\}.*", "\\1", desc)
    probs = strsplit(probs, ",")[[1]]

    if (length(probs) != nlevels)
      stop("the dimension of the CPT of node ", node,
        " does not match the number of its levels.")

    node.cpt = as.table(as.numeric(probs))
    dimnames(node.cpt) = list(nodes.levels[[node]])

  }#THEN
  else {

    row = strsplit(sub(".+\\{[^(]*(\\(.+?)\\s*[;]*\\s*\\}.*", "\\1", desc), ";")[[1]]

    # extract the conditional probability distributions.
    probs = strsplit(sub(".*\\)\\s*[:]{0,1}\\s*(.+)", "\\1", row), ",")
    probs = lapply(probs, as.numeric)
    # paper over degenerate distributions with a uniform distribution.
    probs = lapply(probs, function(x) {

      if (all(x == 0))
        return(prop.table(rep(1, nlevels)))
      else
        return(x)

    })

    # extract the configurations of the parents.
    cfg = strsplit(sub(".*\\((.+)\\).+", "\\1", row), ",")
    cfg = lapply(cfg, sub, pattern = "^\\s*(.+?)\\s*$", replacement = "\\1")
    dims = lapply(c(node, parents[[node]]), function(x) nodes.levels[[x]])

    # DSC files use numeric coordinates instead of levels; Genie does not write
    # them down, it fills them up with zeroes. Use my best guess of the right
    # ordering (most significant corrdinate is the right-most, indexes start
    # from zero) and hope for the best.
    if (all(sapply(cfg, function(x) all(x == 0)))) {

      cfg = lapply(dims, function(x) seq(length(x)))
      cfg = expand.grid(cfg[1+rev(1:length(parents[[node]]))])
      cfg = cfg[rev(1:ncol(cfg))]
      cfg = lapply(seq(nrow(cfg)), function(x) cfg[x, ])

    }#THEN
    else {

      cfg = lapply(cfg, function(x) as.numeric(x) + 1)

    }#ELSE

    # check whether the number of conditional probability distributions matches
    # the number of configurations.
    if (length(probs) != length(cfg))
      stop("the number of conditional distributions for node ", node,
        " do not math the number of configurations of its parents")
    # check whether the number of configurations matches the expected number of
    # configurations given the parents' number of levels.
    expected.cfgs = prod(sapply(nodes.levels[parents[[node]]], length))
    if (length(probs) != expected.cfgs)
      stop("there should be ", expected.cfgs, " conditional probability ",
        "distributions, but only ", length(cfg), " are present.")
    # check whether each conditional probability distribution is valid.
    not.prob = !sapply(probs, is.probability.vector)
    if (any(not.prob)) {

      not.prob = paste(which(not.prob), collapse = ", ")
      stop("conditional probability ditribution(s) ", not.prob, " of node ",
        node, " is not a vector of probabilities.")

    }#THEN
    # check whether each conditional probability distribution sums to one.
    not.prob = (abs(sapply(probs, sum) - 1) > 0.011)
    if (any(not.prob)) {

      # if more than 1% of probability mass, let's assume that the conditional
      # probability distribution is misspecied.
      not.prob = paste(which(not.prob), collapse = ", ")
      stop("conditional probability ditribution(s) ", not.prob, " of node ",
        node, " does not sum to one.")

    }#THEN
    else {

      # perform some more rounding to make the total probability mass closer
      # to one.
      probs = lapply(probs, prop.table)

    }#ELSE
    # check whether the conditional probability distributions have the right
    # number of elements
    if (any(lapply(probs, function(x) length(x)) != nlevels))
      stop("one of the conditional probability ditributions of node ", node,
        " has the wrong number of elements.")

    node.cpt = table(seq(prod(sapply(dims, length))))
    dim(node.cpt) = sapply(dims, length)
    dimnames(node.cpt) = dims

    for (i in seq_along(probs)) {

      node.cpt = do.call("[<-", c(list(node.cpt, 1:nlevels),
                   cfg[[i]], list(probs[[i]])))

    }#FOR

  }#ELSE

  return(node.cpt)

}#DSC.GET.PROBABILITIES

# build the (conditional) probability table for a node in a NET file.
net.get.probabilities = function(node, start, lines, nodes.levels, parents, root) {

  end.line = 0

  # get the line in which the description is starting.
  start.line = start[node]
  # loop until the line in which the description is ending.
  end.line = match.brace(lines, start.line)

  # get all the node's description on one line for easy handling.
  desc = paste(lines[start.line:end.line], collapse = " ")

  # check bogus potential entries
  if (length(grep("model_data\\s+=", desc)) == 1)
    stop("unknown CPT format 'model_data'.")

  # deparse the node's probability table.
  if (root) {

    probs = sub(".+\\{\\s*data\\s*=\\s*\\(\\s*(.+?)\\s*\\);\\s*\\}.*", "\\1", desc)
    probs = gsub("\\s*\\(\\s*|\\s*\\)\\s*", "", probs)
    probs = strsplit(probs, " ")[[1]]

    if (length(probs) != length(nodes.levels[[node]]))
      stop("the dimension of the CPT of node ", node,
        " does not match the number of its levels.")

    node.cpt = as.table(as.numeric(probs))
    dimnames(node.cpt) = list(nodes.levels[[node]])

  }#THEN
  else {

    row = gsub("\\(|\\)", " ", sub(".+\\{[^(]*(\\(.+?)\\s*[;]*\\s*\\}.*", "\\1", desc))

    # extract the conditional probability distributions.
    dims = lapply(c(node, parents[[node]]), function(x) nodes.levels[[x]])
    nconfigs = prod(sapply(dims[-1], length))
    probs = strsplit(sub("^\\s+", "", row), "\\s+")[[1]]

    # check whether we have the expected number of conditional probabilities.
    if (length(probs) != prod(sapply(dims, length)))
      stop("one of the conditional probability ditributions of node ", node,
           " has the wrong number of elements.")

    probs = matrix(as.numeric(probs), byrow = TRUE, nrow = nconfigs)

    # check whether each conditional probability distribution is valid.
    not.prob = !apply(probs, 1, is.probability.vector)
    if (any(not.prob)) {

      not.prob = paste(which(not.prob), collapse = ", ")
      stop("conditional probability ditribution(s) ", not.prob, " of node ",
        node, " is not a vector of probabilities.")

    }#THEN
    # check whether each conditional probability distribution sums to one.
    not.prob = (abs(apply(probs, 1, sum) - 1) > 0.011)
    if (any(not.prob)) {

      # if more than 1% of probability mass, let's assume that the conditional
      # probability distribution is misspecied.
      not.prob = paste(which(not.prob), collapse = ", ")
      stop("conditional probability ditribution(s) ", not.prob, " of node ",
        node, " does not sum to one.")

    }#THEN
    else {

      # perform some more rounding to make the total probability mass closer
      # to one.
      probs = prop.table(probs, margin = 1)

    }#ELSE

    node.cpt = table(seq(prod(sapply(dims, length))))
    dim(node.cpt) = sapply(dims, length)
    dimnames(node.cpt) = dims

    # NET files list conditional probability distributions without mentioning
    # their coordinates in the CPT. Use my best guess of the right ordering
    # (most significant corrdinate is the right-most, indexes start from zero)
    # and hope for the best.
    cfg = lapply(dims, function(x) seq(length(x)))
    cfg = expand.grid(cfg[1+rev(1:length(parents[[node]]))])
    cfg = cfg[rev(1:ncol(cfg))]
    cfg = lapply(seq(nrow(cfg)), function(x) cfg[x, ])

    for (i in seq(nrow(probs))) {

      node.cpt = do.call("[<-", c(list(node.cpt, 1:length(nodes.levels[[node]])),
                   cfg[[i]], list(probs[i, ])))

    }#FOR

  }#ELSE

  return(node.cpt)

}#NET.GET.PROBABILITIES

# match braces to delimit blocks.
match.brace = function(lines, start, open = "{", close = "}") {

  .Call("match_brace",
        lines = lines,
        start = start,
        open = open,
        close = close,
        PACKAGE = "bnlearn")

}#MATCH.BRACE
