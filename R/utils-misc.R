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

  .Call("r_subsets",
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

  .Call("configurations",
        data = data,
        factor = factor,
        all = all)

}#CONFIGURATIONS

# rbind-like function for arc sets.
arcs.rbind = function(matrix1, matrix2, reverse2 = FALSE) {

  .Call("arcs_rbind",
        matrix1 = matrix1,
        matrix2 = matrix2,
        reverse2 = reverse2)

}#ARCS.RBIND

minimal.data.frame = function(lst) {

  .Call("minimal_data_frame",
        obj = lst)

}#MINIMAL.DATA.FRAME

minimal.data.frame.column = function(dataframe, column, drop = TRUE) {

  .Call("dataframe_column",
        dataframe = dataframe,
        column = column,
        drop = drop)

}#MINIMAL.DATA.FRAME.COLUMN

minimal.qr.matrix = function(dataframe, column) {

  .Call("qr_matrix",
        dataframe = dataframe,
        column = column)

}#MINIMAL.QR.MATRIX

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

    if(!any(sapply(l, is.recursive)))
      break
    else
      l = unlist(lapply(l, as.list))

  }#REPEAT

  return(sapply(l, as.character))

}#EXPLODE

# normalize a conditional probability table.
normalize.cpt = function(x) {

  .Call("normalize_cpt",
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
cgsd = function(resid, configs, p) {

  resid = noattr(resid, ok = character(0))
  mult = function(x, p)
           ifelse(missing(p), 1, sqrt((length(x) - 1) / (length(x) - p)))

  if (missing(configs) || is.null(configs)) {

    sd = sd(resid) * mult(resid, p)

  }#THEN
  else {

    sd = by(data = resid, INDICES = configs,
             FUN = function(x) { sd(x) * mult(x, p) })

  }#ELSE

  # replace NA values with zeroes.
  sd[is.na(sd)] = 0

  return(noattr(sd))

}#CGSD

# wrapper around coefficients() to avoid dispatch.
minimal.coefficients = function(x) {

  if (is(x, "penfit"))
    c(x@unpenalized, x@penalized)
  else
    coefficients(x)

}#MINIMAL.COEFFICIENTS

# wrapper around residuals() to avoid dispatch.
minimal.residuals = function(x) {

  if (is(x, "penfit"))
    x@residuals
  else
    residuals(x)

}#MINIMAL.RESIDUALS

# wrapper around fitted() to avoid dispatch.
minimal.fitted = function(x) {

  if (is(x, "penfit"))
    x@fitted
  else
    fitted(x)

}#MINIMAL.FITTED


