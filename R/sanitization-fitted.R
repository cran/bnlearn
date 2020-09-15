
# check the consistency of a bn.fit.dnode or a bn.fit.onode.
check.dnode = function(x, node) {

  # the CPT must be a table.
  if (!is.ndmatrix(x))
    stop("the conditional probability distribution of node ", node,
      " must be a table, a matrix or a multidimensional array.")
  # all elements must be probabilities.
  if (!is.probability.vector(x))
    stop("the probabilities provided for node ",
      node, " do not form a valid conditional probability distribution.")

  # convert the CPT into a table object.
  x = as.table(x)

  # compute the dimensions of the table
  dims = dim(x)
  ndims = length(dims)
  # flatten 1xc tables into 1-dimensional ones.
  if ((ndims == 2) && (dims[1] == 1))
    x = flatten.2d.table(x)
  # update dims and ndims.
  dims = dim(x)
  ndims = length(dims)

  # normalization.
  if (ndims == 1) {

    if (abs(sum(x) - 1) > 0.01)
      stop("the probability distribution of node ", node, " does not sum to one.")

    x = x / sum(x)

  }#THEN
  else {

    check.sum = sapply(margin.table(x,  seq(from = 2, to = ndims)),
                  function(x) abs(sum(x) - 1) > 0.01)

    if (any(check.sum))
      stop("some conditional probability distributions of node ", node, " do not sum to one.")

    x = prop.table(x, margin = seq(ndims)[-1])

  }#ELSE

  return(x)

}#CHECK.DNODE

# check one bn.fit.dnode or bn.fit.onode against another, using only the
# conditional probability table.
check.dnode.vs.dnode = function(new, old) {

  ndims = length(dim(new))

  if (!is(old, c("bn.fit.dnode", "bn.fit.onode")))
    stop("both nodes should be discrete nodes.")

  # now that we are sure that the dimensions are the right ones, reorder
  # them to follow the ordering of the parents in the network.
  ndnn = dimnames(new)
  odnn = dimnames(old$prob)

  if (!is.null(ndnn)) {

    # if the new CPT has dimension names, they must match (and be in the same
    # order too).
    if (ndims > 1) {

      if (!setequal(names(odnn), names(ndnn)))
        stop("wrong dimension names for node ", old$node, ".")

      new = aperm(new, match(names(ndnn), names(odnn)))
      ndnn = dimnames(new)

    }#THEN

    # check labels in the same way, if present.
    if (any(!sapply(ndnn, is.null))) {

        dnn.check = sapply(seq_along(ndnn),
                      function(i) setequal(ndnn[[i]], odnn[[i]]))

        if (!all(dnn.check))
          stop("wrong levels for node ", old$node, ".")

        new = ndsubset(new, odnn)

    }#THEN
    else {

      if (!identical(dim(new), dim(old$prob)))
        stop("wrong dimensions for node ", old$node, ".")

    }#ELSE

  }#THEN
  else {

    # same dimensions (there is nothing else to check).
    if (!identical(dim(new), dim(old$prob)))
      stop("wrong dimensions for node ", old$node, ".")

  }#ELSE

  # make sure the returned object is a table.
  return(structure(new, class = "table"))

}#CHECK.DNODE.VS.DNODE

# check a bn.fit.dnode or a bn.fit.onode against the local distributions of
# its parents.
check.dnode.vs.parents = function(node, new, parents) {

  if (length(parents) == 0)
    return(new)

  # all parents of discrete nodes must be discrete nodes themselves.
  illegal.parents =
    sapply(parents, class) %!in% c("bn.fit.dnode", "bn.fit.onode")

  if (any(illegal.parents))
    stop("node ", node, " is discrete but has continuous parents (",
      names(parents)[illegal.parents], ").")
  # extract the CPTs.
  cpt.levels = dimnames(new)
  parents.levels = smartSapply(NULL, names(parents),
                     function(p) dimnames(parents[[p]]$prob)[[1]])

  # check the number of dimensions.
  if (length(cpt.levels) != length(parents) + 1)
    stop("wrong number of parents for node ", node, ".")
  # check the dimension names, and add them if they are missing
  if (is.null(names(cpt.levels)))
    names(dimnames(new)) = names(cpt.levels) = c(node, names(parents.levels))
  else {

    # check whether the CPT covers all the relevant variables.
    if (!setequal(names(cpt.levels)[-1], names(parents.levels)))
      stop("wrong dimensions for node ", node, ".")
    # now that we are sure that the dimensions are the right ones, reorder
    # them to follow the ordering of the parents in the network.
    d = names(dimnames(new))
    new = aperm(new, c(match(node, d), match(names(parents.levels), d)))
    cpt.levels = dimnames(new)

  }#ELSE

  # check the levels in the CPT against those of the parents, and reoder them
  # if needed.
  for (p in names(parents))
    if (!setequal(parents.levels[[p]], cpt.levels[[p]]))
      stop("the levels of the parent ", p, " of node ", node, " do not match.")

  if (!isTRUE(all.equal(parents.levels, cpt.levels[-1])))
    new = ndsubset(new, parents.levels)

  return(new)

}#CHECK.DNODE.VS.PARENTS

# check the consistency of a bn.fit.gnode or a bn.fit.cgnode.
check.gnode = function(x, node) {

  components =  c("coef", "fitted", "resid", "sd", "configs")
  labels = c(fitted = "fitted values", resid = "residuals")

  # custom list of components.
  if (!is.list(x) || any(names(x) %!in% components))
    stop("the conditional probability distribution for node ", node,
         " must be a list with at least one of the following elements:",
         paste0(" '", components, "'"), ".")
  if (("coef" %!in% names(x)) || !any(c("sd", "resid") %in% names(x)))
    stop("at least the regression coefficients and either the residuals or ",
      "the residual standard deviation are required for node ", node, ".")
  # discrete parents' configurations only belong to bn.fit.cgnode.
  if ((length(dim(x$coef)) != 2) && ("configs" %in% names(x)))
    stop("no parents' configurations are needed with a single set of coefficients.")

  if (!is.null(x$coef))
    if ((length(x$coef) == 0) || !is.real.vector(x$coef))
      stop("coef must be a vector or a matrix of numeric values, the ",
        "regression coefficients for node ", node, " given its parents.")

  for (comp in c("fitted", "resid"))
    if (!is.null(x[[comp]]))
      if ((length(x[[comp]]) == 0) || !is.real.vector(x[[comp]]))
        stop(comp, " must be a vector of numeric values, the ",
          labels[comp], " for node ", node, " given its parents.")

  if (!is.null(x$configs)) {

    if ((length(x$configs) == 0) || !is.factor(x$configs))
      stop("the discrete parents' configurations for node ", node,
        " must be a factor.")
    if (nlevels(x$config) != ncol(x$coef))
      stop("wrong number of discrete parents' configurations for node ", node, ".")

  }#THEN

  # check the standard deviation of the residuals.
  if (!is.null(x$sd)) {

    if (!is.nonnegative.vector(x$sd))
      stop("sd must be non-negative, the standard deviation(s) of the ",
        "residuals of node ", node, ".")

    if ((length(dim(x$coef)) == 2) && (length(x$sd) != ncol(x$coef)) ||
        (length(dim(x$sd)) > 1))
      stop("the dimensions of sd and coef do not match.")
    if ((length(dim(x$coef)) %in% c(0, 1)) && (length(x$sd) > 1))
      stop("sd chould be a single non-negative number.")

    if (!is.null(x$resid) && (length(x$sd) == 1)) {

      adj.sd = cgsd(x$resid, p = length(x$coef))

      if (!isTRUE(all.equal(x$sd, adj.sd, check.attributes = FALSE, tol = 0.0005)))
        stop("the reported standard deviation of the residuals of node ", node,
          " does not match the observed one.")

    }#THEN
    else if (!is.null(x$resid) && (length(x$sd) > 1) && !is.null(x$configs)) {

      adj.sd = cgsd(x$resid, configs = x$configs, p = nrow(x$coef))

      if (!isTRUE(all.equal(x$sd, adj.sd, check.attributes = FALSE, tol = 0.0005)))
        stop("the reported standard deviation of the residuals of node ", node,
          " does not match the observed one.")

    }#THEN

  }#THEN
  else if (!is.null(x$resid)) {

    # compute the standard error from the residuals if possible.
    if ((length(dim(x$coef)) == 2) && is.null(x$config))
      stop("sd is missing, and parents' configurations are required to compute it.")

    x$sd = cgsd(x$resid, configs = x$configs,
             p = ifelse(is.matrix(x$coef), nrow(x$coef), length(x$coef)))

  }#THEN

  # one residual for each fitted value.
  if (!is.null(x$resid) && !is.null(x$fitted))
    if (length(x$resid) != length(x$fitted))
      stop("the residuals and the fitted values of node ", node,
        " have different lengths.")
  # if any, one parent configuration for each fitted value and residual.
  if (!is.null(x$config)) {

    if (!is.null(x$resid) && (length(x$configs) != length(x$resid)))
      stop("parents' configurations and residuals of node ", node,
        " have different lengths.")
    if (!is.null(x$fitted) && (length(x$configs) != length(x$fitted)))
      stop("parents' configurations and fitted values of node ", node,
        " have different lengths.")

  }#THEN

  return(x)

}#CHECK.GNODE

# check one bn.fit.gnode against another.
check.gnode.vs.gnode = function(new, old) {

  if (!is(old, "bn.fit.gnode"))
    stop("both nodes should be continuous (Gaussian) nodes.")

  # same number of coefficients.
  if (length(new$coef) != length(old$coefficients))
    stop("wrong number of coefficients for node ", old$node, ".")
  # same regression coefficients.
  if (!is.null(names(new$coef))) {

    # if the new coefficients have labels, they must match.
    check.nodes(names(new$coef), graph = names(old$coefficients),
      min.nodes = length(names(old$coefficients)))
    # reorder the coefficients to match.
    new$coef = new$coef[names(old$coefficients)]

  }#THEN
  else {

    # copy the names from the spec.
    names(new$coef) = names(old$coefficients)

  }#ELSE
  # same number of standard errors (fixed to 1).
  if (length(new$sd) != 1)
    stop("wrong number of standard errors for node ", old$node, ".")
  # same number of residuals.
  if (!is.null(new$resid))
    if (length(new$resid) != length(old$residuals))
      stop("wrong number of residuals for node ", old$node, ".")
  # same number of fitted values.
  if (!is.null(new$fitted))
    if (length(new$fitted) != length(old$fitted.values))
      stop("wrong number of fitted values for node ", old$node, ".")

  return(new)

}#CHECK.GNODE.VS.GNODE

# check a bn.fit.gnode using only the names of its parents.
check.gnode.vs.parents = function(node, new, parents) {

  # all parents of Gaussian nodes must be continuous nodes.
  illegal.parents =
    sapply(parents, class) %!in% c("bn.fit.gnode", "bn.fit.cgnode")

  if (any(illegal.parents))
    stop("node ", node, " is Gaussian but has discrete parents (",
      names(parents)[illegal.parents], ").")

  # add the intercept, which is obviously not among the parents of the node.
  parent.names = c("(Intercept)", names(parents))

  # same number of coefficients.
  if (length(new$coef) != length(parent.names))
    stop("wrong number of coefficients for node ", node, ".")
  # same number of standard errors (fixed to 1).
  if (length(new$sd) != 1)
    stop("wrong number of standard errors for node ", node, ".")
  # if the new coefficients have labels, they must match.
  if (!is.null(names(new$coef))) {

    if (!setequal(names(new$coef), parent.names))
      stop("wrong regression coefficients for node ", node, " (",
        paste(setdiff(parent.names, names(new$coef)), collapse = " "), ").")
    new$coef = new$coef[parent.names]

  }#THEN
  else {

    names(new$coef) = parent.names

  }#ELSE

  return(new)

}#CHECK.GNODE.VS.PARENTS

# check one bn.fit.cgnode against another.
check.cgnode.vs.cgnode = function(new, old) {

  if (!is(old, "bn.fit.cgnode"))
    stop("both nodes should be continuous (conditional Gaussian) nodes.")

  # right dimensions for the coefficients.
  if (!is(new$coef, "matrix") ||
      !identical(dim(new$coef), dim(old$coefficients)))
    stop("the regression coefficients for node ", old$node, " must be ",
         "in a matrix with one column for each discrete parent and one ",
         "coefficient for each continuous parent.")
  # right row names for the coefficients.
  if (!is.null(rownames(new$coef))) {

    # if the new coefficients have labels, they must match.
    check.nodes(rownames(new$coef), graph = rownames(old$coefficients),
      min.nodes = length(rownames(old$coefficients)))
    # reorder the coefficients to match.
    new$coef = new$coef[rownames(old$coefficients), , drop = FALSE]

  }#THEN
  else {

    # copy the names from the spec.
    rownames(new$coef) = rownames(old$coefficients)

  }#ELSE
  # right column names for the coefficients.
  if (!is.null(colnames(new$coef))) {

    # if the new coefficients have labels, they must match.
    check.nodes(colnames(new$coef), graph = colnames(old$coefficients),
      min.nodes = length(colnames(old$coefficients)))
    # reorder the coefficients to match.
    new$coef = new$coef[, colnames(old$coefficients), drop = FALSE]

  }#THEN
  else {

    # copy the names from the spec.
    colnames(new$coef) = colnames(old$coefficients)

  }#ELSE
  # same number of residuals.
  if (!is.null(new$resid))
    if (length(new$resid) != length(old$residuals))
      stop("wrong number of residuals for node ", old$node, ".")
  # same number of fitted values.
  if (!is.null(new$fitted))
    if (length(new$fitted) != length(old$fitted.values))
      stop("wrong number of fitted values for node ", old$node, ".")

  if (!is.null(new$dlevels)) {

    old.dparents = names(old$dlevels)
    new.dparents = names(new$dlevels)

    # different sets of discrete parents.
    if (!setequal(old.dparents, new.dparents))
      stop("different discrete parents for node ", old$node, ".")
    # identical sets of discrete parents, but in a different order; match
    # configurations.
    old.config = old$dlevels
    new.config = new$dlevels
    for (i in seq_along(old.config))
      old.config[[i]] = paste(old.dparents[i], old.config[[i]], sep = ":")
    for (i in seq_along(new.config))
      new.config[[i]] = paste(new.dparents[i], new.config[[i]], sep = ":")

    old.config = interaction(expand.grid(old.config))
    new.config = interaction(expand.grid(new.config)[old.dparents])
    map = match(old.config, new.config)

    # reorder coefficients, standard errors and configurations.
    new$sd = new$sd[map]
    new$coef = new$coef[, map, drop = FALSE]
    new$dlevels = new$dlevels[old.dparents]
    colnames(new$coef) = names(new$sd) = names(old$sd)

  }#THEN

  return(new)

}#CHECK.CGNODE.VS.CGNODE

# check a bn.fit.cgnode against the local distributions of its parents.
check.cgnode.vs.parents = function(node, new, parents) {

  # both continuous and discrete parents are legal, but treated differently.
  discrete.id = sapply(parents, class) %in% c("bn.fit.dnode", "bn.fit.onode")
  discrete.parents = names(parents)[discrete.id]
  continuous.id = sapply(parents, class) %in% c("bn.fit.gnode", "bn.fit.cgnode")
  continuous.parents = names(parents)[continuous.id]
  parent.configurations =
    prod(sapply(parents[discrete.parents], function(p) dim(p$prob)[1]))
  correct.configs = as.character(seq(from = 0, to = parent.configurations - 1))

  new$dparents = which(discrete.id)
  new$gparents = which(continuous.id)
  new$dlevels = lapply(parents[discrete.parents], function(p) dimnames(p$prob)[[1]])

  # the node should have at least one discrete parent, or it would not be CLG.
  if (length(discrete.parents) == 0)
    stop("node ", node, " is a conditioanl Gaussian node but has no discrete parent.")

  # check the dimensions of the regression coefficients matrix.
  if (!isTRUE(all.equal(dim(new$coef), c(length(continuous.parents) + 1L,
      parent.configurations))))
    stop("wrong dimensions for the regression coefficients matrix of node ",
      node, ".")

  # check and fix the coefficients names.
  if (!is.null(rownames(new$coef))) {

    regression.coefficients = c("(Intercept)", continuous.parents)

    if (!setequal(rownames(new$coef), regression.coefficients))
      stop("wrong regression coefficients for node ", node, " (",
        paste(setdiff(continuous.parents, rownames(new$coef)), collapse = " "), ").")

    new$coef = new$coef[regression.coefficients, , drop = FALSE]

  }#THEN
  else {

    rownames(new$coef) = c("(Intercept)", continuous.parents)

  }#ELSE

  # check and fix the column names of the coefficients matrix.
  if (is.null(colnames(new$coef)))
    colnames(new$coef) = correct.configs
  else if (setequal(colnames(new$coef), correct.configs))
    new$coef = new$coef[, correct.configs, drop = FALSE]
  else {

    # wrong names: warn and rename.
    warning("remapping levels of the discrete parents configurations ",
      "for node ", node, ".")

    colnames(new$coef) = correct.configs

  }#ELSE

  # check the number of standard errors, and check their names.
  if (length(new$sd) != parent.configurations)
    stop("wrong number of standard errors for node ", node, ".")
  if (is.null(names(new$sd)))
    names(new$sd) = correct.configs
  else if (setequal(names(new$sd), correct.configs))
    new$sd = new$sd[correct.configs]
  else {

    # wrong names: warn and rename.
    warning("remapping levels of the discrete parents configurations ",
      "for node ", node, ".")

    names(new$sd) = correct.configs

  }#ELSE

  # check the discrete parents configurations.
  if (!is.null(new$configs)) {

    actual.configs = levels(new$configs)

    if (length(actual.configs) != parent.configurations)
      stop("wrong number of discrete parents' configurations for node", node, ".")
    if (any(actual.configs != correct.configs)) {

      if (setequal(actual.configs, correct.configs))
        new$configs = factor(new$configs, correct.configs)
      else {

        # wrong levels: warn and rename.
        warning("remapping levels of the discrete parents configurations ",
          "for node ", node, ".")

        levels(new$configs) = correct.configs

      }#ELSE

    }#THEN

  }#THEN

  return(new)

}#CHECK.CGNODE.VS.PARENTS

# check discrete parents configurations in list form for all kinds of nodes.
check.discrete.parents.configuration = function(config, node,
    ideal.only = FALSE) {

  # check whether the configuration is there.
  if (missing(config) || !is(config, "list"))
    stop("the parents configuration must be a list with elements named after the parents of ", node$node, ".")

  if (is(node, c("bn.fit.dnode", "bn.fit.onode")))
    allowed.parents = node$parents
  else if (is(node, c("bn.fit.cgnode")))
    allowed.parents = node$parents[node$dparents]

  # check the node labels from the list names.
  check.nodes(names(config), graph = allowed.parents,
    min.nodes = length(allowed.parents), max.nodes = length(allowed.parents))

  # check that only one value is provided for each parent, to get a single
  # configuration.
  for (c in names(config))
    if (length(config[[c]]) > 1)
      stop("only one value allowed for parent '", c, "'.")

  if (is(node, c("bn.fit.dnode", "bn.fit.onode")))
    levels = dimnames(node$prob)
  else
    levels = node$dlevels

  # check whether the provided values are valid values.
  for (c in names(config)) {

    if (is.factor(config[[c]]))
      config[[c]] = as.character(config[[c]])

    if (config[[c]] %!in% levels[[c]])
      stop("the value of parent '", c, "' must be a valid level.")

  }#FOR

  return(config)

}#CHECK.DISCRETE.PARENTS.CONFIGURATION

