
# check the data set.
check.data = function(x, allowed.types = available.data.types,
    allow.levels = FALSE, allow.missing = FALSE, warn.if.no.missing = FALSE,
    stop.if.all.missing = FALSE, stop.if.missing.whole.obs = FALSE) {

  # check the data are there.
  if (missing(x))
    stop("the data are missing.")
  # x must be a data frame.
  if (!is.data.frame(x))
    stop("the data must be in a data frame.")
  # x must be a data frame with at least one column.
  if (ncol(x) == 0)
    stop("the data must be in a data frame with at least one column.")
  # check the data for NULL/NaN/NA.
  observed = count.observed.values(x)

  if (allow.missing) {

   if (warn.if.no.missing && all(observed$columns == nrow(x)))
     warning("no missing data are present even though some are expected.")

   if (any(observed$columns == 0))
     if (stop.if.all.missing)
       stop("at least one variable has no observed values.")
     else
       warning("at least one variable has no observed values.")

  }#THEN
  else {

    if (any(observed$columns < nrow(x)))
      stop("the data set contains NaN/NA values.")

  }#ELSE

  if (any(observed$rows == 0))
    if (stop.if.missing.whole.obs)
      stop("some observations contain only missing values.")
    else
      warning("some observations contain only missing values.")

  # check which type of data we are dealing with.
  type = data.type(x)
  # check whether the variables are of the expected types.
  check.label(type, choices = allowed.types, labels = data.type.labels,
    argname = "data type")

  # checks specific to a particular data type.
  if (type %in% discrete.data.types) {

    for (col in names(x)) {

      # check the number of levels of discrete variables, to guarantee that
      # the degrees of freedom of the tests are positive.
      if (nlevels(x[, col]) < 2)
        stop("variable ", col, " must have at least two levels.")

      # warn about levels with zero frequencies, it's not necessarily wrong
      # (data frame subsetting) but sure is fishy.
      counts = .table(x[, col, drop = FALSE], with.missing = allow.missing)
      if (!allow.levels && any(counts == 0))
        warning("variable ", col, " has levels that are not observed in the data.")

    }#FOR

  }#THEN
  else if (type == "continuous") {

    # all values must be finite to have finite mean and variance.
    check.data.frame.finite(x)

  }#THEN

  return(list(type = type, complete.nodes = (observed$columns == nrow(x))))

}#CHECK.DATA

# is the data of a particular type?
data.type = function(data) {

  .Call(call_data_type,
        data = data)

}#DATA.TYPE

# count observed values per-variable and per-observation in a single pass.
count.observed.values = function(data) {

  .Call(call_count_observed_values,
        data = data)

}#COUNT.OBSERVED.VALUES
