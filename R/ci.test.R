
# do a single conditional independence test.
ci.test = function(x, ...) {

  UseMethod("ci.test", x)

}#CI.TEST

# do a single conditional independence test (nodes as character strings).
ci.test.character = function(x, y = NULL, z = NULL, data, test = NULL,
    B = NULL, debug = FALSE, ...) {

  # the original data set is needed.
  check.data(data)
  # check debug.
  check.logical(debug)
  # check the variables involved in the test.
  if (!is.string(x) || !(x %in% names(data)))
    stop("'x' must be a character string, the name of one of the columns of 'data'.")
  if (!is.string(y) || !(y %in% names(data)))
    stop("'y' must be a character string, the name of one of the columns of 'data'.")
  if (x == y)
    stop("'x' must be different from 'y'.")
  if (!is.null(z)) {

    if (!is.character(z) || !all(z %in% names(data)))
      stop("'z' must be a vector of character strings, the names of one or more of the columns of 'data'.")
    if (any(z %in% c(x, y)))
      stop("'z' must be different from both 'x' and 'y'.")

  }#THEN
  # check the test label.
  test = check.test(test, data)
  # check B (the number of bootstrap/permutation samples).
  B = check.B(B, test)
  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  # compute the network score.
  conditional.test(x = x, y = y, sx = z, data = data, test = test,
    B = B, alpha = 1, learning = FALSE)

}#CI.TEST.CHARACTER

# do a single conditional independence test (nodes in a matrix).
ci.test.matrix = function(x, test = NULL, B = NULL, debug = FALSE, ...) {

  if (ncol(x) < 2)
    stop("'x' must have at least two columns.")

  ci.test(x = as.data.frame(x), test = test, B = B, debug = debug, ...)

}#CI.TEST.MATRIX

# do a single conditional independence test (nodes in a data frame).
ci.test.data.frame = function(x, test = NULL, B = NULL, debug = FALSE, ...) {

  if (ncol(x) < 2)
    stop("'x' must have at least two columns.")

  nodes = names(x)

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  ci.test.character(x = nodes[1], y = nodes[2], z = nodes[-(1:2)],
    data = x, test = test, B = B, debug = debug)

}#CI.TEST.DATA.FRAME

# do a single conditional independence test (numerical vectors).
ci.test.numeric = function(x, y = NULL, z = NULL, test = NULL, B = NULL, debug = FALSE, ...) {

  # check debug.
  check.logical(debug)
  # check the variables involved in the test.
  if (!is.numeric(y) && !(is.matrix(y) && ncol(y) == 1))
    stop("'y' must be a numeric vector.")
  if (length(y) != length(x))
    stop("'x' and 'y' must have the same length.")
  if (!is.null(z)) {

    if (is.matrix(z) || is.data.frame(z)) {

      if (nrow(z) != length(x))
        stop("'x', 'y', and 'z' must have the same length.")
      if (!is.data.continuous(z))
        stop("'z' must be a numeric matrix or data frame.")

      sx = 3:(2 + ncol(z))

    }#THEN
    else if (is.numeric(z)) {

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

  }#ELSE
  # check the data are there.
  check.data(data)
  # check the test label.
  test = check.test(test, data)
  # check B (the number of bootstrap/permutation samples).
  B = check.B(B, test)
  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  res = conditional.test(x = 1L, y = 2L, sx = sx, data = data,
    test = test, B = B, alpha = 1, learning = FALSE)

  # rewrite the test formula.
  res$data.name = paste(deparse(substitute(x)), "~", deparse(substitute(y)),
        ifelse(!is.null(z), paste("|", deparse(substitute(z))), ""),
        collapse = " + ")

  return(res)

}#CI.TEST.NUMERIC

# do a single conditional independence test (factor objects).
ci.test.factor = function(x, y = NULL, z = NULL, test = NULL, B = NULL,
    debug = FALSE, ...) {

  # check debug.
  check.logical(debug)
  # check the variables involved in the test.
  if (!is.factor(y) && !(is.matrix(y) && ncol(y) == 1))
    stop("'y' must be a factor vector.")
  if (length(y) != length(x))
    stop("'x' and 'y' must have the same length.")
  if (!is.null(z)) {

    if (is.matrix(z) || is.data.frame(z)) {

      if (nrow(z) != length(x))
        stop("'x', 'y', and 'z' must have the same length.")
      if (!is.data.discrete(z))
        stop("'z' must be a factor matrix or data frame.")

      sx = 3:(2 + ncol(z))

    }#THEN
    else if (is.factor(z)) {

      if (length(z) != length(x))
        stop("'x', 'y', and 'z' must have the same length.")

      sx = 3L

    }#THEN
    else
      stop("'z' must be a factor vector, matrix or data frame.")

    # build the data frame.
    data = data.frame(x = x, y = y, z = z)

  }#THEN
  else {

    # build the data frame.
    data = data.frame(x = x, y = y)

    sx = integer(0)

  }#ELSE
  # check the data are there.
  check.data(data)
  # check the test label.
  test = check.test(test, data)
  # check B (the number of bootstrap/permutation samples).
  B = check.B(B, test)
  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  res = conditional.test(x = 1L, y = 2L, sx = as.integer(sx), data = data,
    test = test, B = B, alpha = 1, learning = FALSE)

  # rewrite the test formula.
  res$data.name = paste(deparse(substitute(x)), "~", deparse(substitute(y)),
        ifelse(!is.null(z), paste("|", deparse(substitute(z))), ""),
        collapse = " + ")

  return(res)

}#CI.TEST.FACTOR

ci.test.default = function(x, ...) {

  stop("x must be either a factor object, a numeric vector, a character string or a data frame.")

}#CI.TEST.DEFAULT

