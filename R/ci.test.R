
# do a single conditional independence test.
ci.test = function(x, y, z, data, test, B, debug = FALSE) {

  if (missing(x))
    stop("one or both of the variables to test are missing.")

  if (is.string(x)) {

    ci.test.character(x = x, y = y, z = z, data = data, test = test, B = B,
      debug = debug)

  }#THEN
  else if (is(x, c("matrix", "data.frame"))) {

    if (!missing(y) || !missing(z))
      warning("'y' and 'z' will be ignored.")
    if (!missing(data))
      warning("'data' will be ignored.")
    if (ncol(x) < 2)
      stop("'x' must have at least two columns.")

    nodes = names(x)
    ci.test.character(x = nodes[1], y = nodes[2], z = nodes[-(1:2)],
      data = x, test = test, B = B, debug = debug)

  }#THEN
  else if (is.vector(x)) {

    if (!missing(data))
      warning("'data' will be ignored.")

    ci.test.vector(x = x, y = y, z = z, xlab = deparse(substitute(x)),
      ylab = deparse(substitute(y)), zlab = deparse(substitute(z)),
      test = test, B = B, debug = debug)

  }#THEN
  else {

    stop("x must be either a factor object, a numeric vector, a character string or a data frame.")

  }#ELSE

}#CI.TEST

# do a single conditional independence test (nodes as character strings).
ci.test.character = function(x, y, z, data, test, B, debug = FALSE) {

  # the original data set is needed.
  check.data(data)
  # check debug.
  check.logical(debug)
  # check the variables involved in the test.
  if (missing(x) || missing(y))
    stop("one or both of the variables to test are missing.")
  if (!is.string(x) || (x %!in% names(data)))
    stop("'x' must be a character string, the name of one of the columns of 'data'.")
  if (!is.string(y) || (y %!in% names(data)))
    stop("'y' must be a character string, the name of one of the columns of 'data'.")
  if (x == y)
    stop("'x' must be different from 'y'.")
  if (!missing(z) && !identical(z, character(0))) {

    if (!is.string.vector(z) || any(z %!in% names(data)))
      stop("'z' must be a vector of character strings, the names of one or more of the columns of 'data'.")
    if ((z == x) || (z == y))
      stop("'z' must be different from both 'x' and 'y'.")

  }#THEN
  else {

    z = character(0)

  }#ELSE
  # check the test label.
  test = check.test(test, minimal.data.frame.column(data, c(x, y, z)))
  # check B (the number of permutation samples).
  B = check.B(B, test)

  # create the htest object.
  htest = indep.test(x = x, y = y, sx = z, data = data, test = test,
            B = B, alpha = 1, learning = FALSE)
  htest$method = test.labels[test]
  htest$data.name = paste(x, "~", y, ifelse(length(z) > 0, "|", ""),
                      paste(z, collapse = " + "))

  return(htest)

}#CI.TEST.CHARACTER

# do a single conditional independence test (data vectors).
ci.test.vector = function(x, y, z, xlab, ylab, zlab, test, B, debug = FALSE) {

  # check debug.
  check.logical(debug)
  # check the variables involved in the test.
  if (!is.vector(x))
    stop("'y' must be a vector.")
  if (!is.vector(y))
    stop("'y' must be a vector.")
  if (length(y) != length(x))
    stop("'x' and 'y' must have the same length.")
  if (!missing(z)) {

    if (is.matrix(z) || is.data.frame(z)) {

      if (nrow(z) != length(x))
        stop("'x', 'y', and 'z' must have the same length.")

      sx = 3:(2 + ncol(z))

    }#THEN
    else if (is.vector(z)) {

      if (length(z) != length(x))
        stop("'x', 'y', and 'z' must have the same length.")

      sx = 3L

    }#THEN
    else
      stop("'z' must be a numeric vector, matrix or data frame.")

    # build the data frame.
    data = data.frame(x = x, y = y, z = z)

  }#THEN
  else {

    # build the data frame.
    data = data.frame(x = x, y = y)
    sx = integer(0)
    z = character(0)

  }#ELSE

  # check the data are there.
  check.data(data)
  # check the test label.
  test = check.test(test, data)
  # check B (the number of permutation samples).
  B = check.B(B, test)

  # create the htest object.
  htest = indep.test(x = 1L, y = 2L, sx = sx, data = data,
    test = test, B = B, alpha = 1, learning = FALSE)
  htest$method = test.labels[test]
  htest$data.name = paste(xlab, "~", ylab,
        ifelse(length(z) > 0, paste("|", zlab), ""))

  return(htest)

}#CI.TEST.VECTOR

