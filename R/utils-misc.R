# this is to keep the old S3 behaviour inside the NAMESPACE.
is = function(x, class) {

  if (identical(class, "double"))
    is.double(x)
  else if ("double" %in% class)
    any(class(x) %in% class) || is.double(x)
  else
    any(class(x) %in% class)

}#IS

# get all the subsets of a given size, even if either the initial set
# or the subset are empty (i.e. of size zero).
subsets = function(elems, size) {

  # allow empty subsets (i.e. subsets of empty sets).
  if ((length(elems) == 0) || (size == 0))
    return(matrix(character(0), nrow = 0, ncol = 0))

  .Call(call_subsets,
        elems = elems,
        size = as.integer(size))

}#SUBSETS

# return the array whose size is smaller.
smaller = function(a, b) {

  if (length(a) < length(b))
    return(a)
  else
    return(b)

}#SMALLER

# build an array containing the configurations of the variables.
configurations = function(data, factor = TRUE, all = TRUE) {

  .Call(call_configurations,
        data = data,
        factor = factor,
        all = all)

}#CONFIGURATIONS

# rbind-like function for arc sets.
arcs.rbind = function(matrix1, matrix2, reverse2 = FALSE) {

  .Call(call_arcs_rbind,
        matrix1 = matrix1,
        matrix2 = matrix2,
        reverse2 = reverse2)

}#ARCS.RBIND

.table = function(x, with.missing = FALSE) {

  .Call(call_minimal_table,
        x = x,
        missing = with.missing);

}#.TABLE

.data.frame = function(lst) {

  .Call(call_minimal_data_frame,
        obj = lst)

}#.DATA.FRAME

.data.frame.column = function(dataframe, column, drop = TRUE) {

  .Call(call_dataframe_column,
        dataframe = dataframe,
        column = column,
        drop = drop)

}#.DATA.FRAME.COLUMN

# flatten 2-dimensional 1xc tables.
flatten.2d.table = function(x) {

  x = as.table(structure(c(x), names = colnames(x)))
  names(dimnames(x)) = ""

  return(x)

}#FLATTEN.2D.TABLE

# explode an unevaluated expression into a character vector.
explode = function(x) {

  if (is.list(x))
    return(names(x))
  else if (identical(x, TRUE))
    return(character(0))
  else
    l = as.list(x)

  repeat {

    if (!any(sapply(l, is.recursive)))
      break
    else
      l = unlist(lapply(l, as.list))

  }#REPEAT

  return(sapply(l, as.character))

}#EXPLODE

# normalize a conditional probability table.
normalize.cpt = function(x) {

  .Call(call_normalize_cpt,
        cpt = x)

}#NORMALIZE.CPT

# negated inclusion operator.
`%!in%` = function(x, table) {

  match(x, table, nomatch = 0L) == 0L

}#%!IN%

# remove extraneous attributes.
noattr = function(x, ok) {

  if (missing(ok))
    if (is.matrix(x))
      ok = c("dim", "dimnames")
    else if (is.factor(x))
      ok = c("class", "levels")
    else
      ok = character(0)

  x.attr = attributes(x)
  attributes(x) = x.attr[names(x.attr) %in% ok]

  return(x)

}#NOATTR

# reset the attributes of a CPT (with dimnames is complicated)
cptattr = function(cpt) {

  # marginal tables have no dimension names (and a single dimension).
  if (length(dim(cpt)) == 1)
    dnn = noattr(dimnames(cpt))
  else
    dnn = dimnames(cpt)
  dim(cpt) = noattr(dim(cpt))
  dimnames(cpt) = dnn

  return(cpt)

}#CPTATTR

# conditional standard deviation.
cgsd = function(x, configs = NULL, p = 1L) {

  .Call(call_cgsd,
        x = x,
        strata = configs,
        nparams = p)

}#CGSD

# wrapper around coefficients() to avoid dispatch.
.coefficients = function(x) {

  if (is(x, "penfit")) {

    # discard zero coefficients from LASSO models.
    if ((x@lambda2 == 0) &&(x@lambda1 > 0))
      c(x@unpenalized, x@penalized[x@penalized != 0])
    else
      c(x@unpenalized, x@penalized)

  }#THEN
  else {

    coefficients(x)

  }#ELSE

}#.COEFFICIENTS

# wrapper around residuals() to avoid dispatch.
.residuals = function(x) {

  if (is(x, "penfit"))
    x@residuals
  else
    residuals(x)

}#.RESIDUALS

# wrapper around fitted() to avoid dispatch.
.fitted = function(x) {

  if (is(x, "penfit"))
    x@fitted
  else
    fitted(x)

}#.FITTED

# subset an n-dimensional matrix in a programmatic way.
ndsubset = function(x, indices) {

  if (length(dim(x)) > 1) {

    index = as.list(structure(rep(TRUE, length(dim(x))), names = names(dimnames(x))))
    index[names(indices)] = indices

  }#THEN
  else {

    index = indices

  }#ELSE

  do.call(`[`, c(list(x), index))

}#NDSUBSET

# make sure rounded probabilites sum up to one (largest remainder method).
lrm.round = function(prob, digits = 3) {

  # scale the probabilities so that the last significant digit is at 10^0.
  scaled = prob * 10^digits
  # separate integer and fractional parts, and the lost probability mass
  integer.part = floor(scaled)
  fractional.part = scaled - integer.part
  # if the resulting number are round, the probabilities are alreay rounded.
  if (isTRUE(all.equal(scaled, integer.part)))
    return(prob)
  # compute the lost probability mass, and where to add it back.
  lost.prob = sum(scaled) - sum(integer.part)
  add.back = order(fractional.part, decreasing = TRUE)[1:lost.prob]
  # add it back.
  integer.part[add.back] = integer.part[add.back] + 1

  # rescale back before returning.
  return(integer.part / 10^digits)

}#LRM.ROUND

