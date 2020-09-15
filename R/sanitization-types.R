
# is x a vector (as opposed to a matrix)?
is.vector = function(x) {

  is.null(dim(x)) ||
  (length(dim(x)) == 2) &&
  (dim(x)[2] == 1)

}#IS.VECTOR

# is x a real number?
is.real.number = function(x) {

  is.numeric(x) &&
  (length(x) == 1) &&
  is.finite(x)

}#IS.REAL

# is x a vector of real number?
is.real.vector = function(x) {

  is.numeric(x) &&
  all(is.finite(x))

}#IS.REAL.VECTOR

# is x a positive number?
is.positive = function(x) {

  is.numeric(x) &&
  (length(x) == 1) &&
  is.finite(x) &&
  (x > 0)

}#IS.POSITIVE

# is x a non-negative number?
is.non.negative = function(x) {

  is.numeric(x) &&
  (length(x) == 1) &&
  is.finite(x) &&
  (x >= 0)

}#IS.NON.NEGATIVE

# is x a positive integer?
is.positive.integer = function(x) {

  is.positive(x) && ((x %/% 1) == x)

}#IS.POSITIVE.INTEGER

# is x a non-negative integer?
is.non.negative.integer = function(x) {

  is.non.negative(x) && ((x %/% 1) == x)

}#IS.NON.NEGATIVE.INTEGER

# is x a vector of positive numbers?
is.positive.vector = function(x) {

  is.numeric(x) &&
  all(is.finite(x)) &&
  all(x > 0)

}#IS.POSITIVE.VECTOR

# is x a vector of non-negative numbers?
is.nonnegative.vector = function(x) {

  is.numeric(x) &&
  all(is.finite(x)) &&
  all(x >= 0)

}#IS.NONNEGATIVE.VECTOR

# is x a probability?
is.probability = function(x) {

  is.numeric(x) &&
  (length(x) == 1) &&
  is.finite(x) &&
  (x >= 0) &&
  (x <= 1)

}#IS.PROBABILITY

# is x a vector of probabilities?
is.probability.vector = function(x, zero = FALSE) {

  is.numeric(x) &&
  all(is.finite(x)) &&
  all(x >= 0) &&
  all(x <= 1) &&
  (zero || any(x > 0))

}#IS.PROBABILITY.VECTOR

# is x a single character string?
is.string = function(x) {

  is.character(x) &&
  (length(x) == 1) &&
  !any(is.na(x)) &&
  any(x != "")

}#IS.STRING

# is x a vector of character strings?
is.string.vector = function(x) {

  is.character(x) &&
  !any(is.na(x)) &&
  any(x != "")

}#IS.STRING.VECTOR

is.logical.vector = function(x) {

  is.logical(x) &&
  !any(is.na(x))

}#IS.LOGICAL.VECTOR

is.ndmatrix = function(x) {

  is(x, c("table", "matrix", "array"))

}#IS.NDMATRIX

# are all numeric values in the data frame fimite?
check.data.frame.finite = function(x) {

  .Call(call_data_frame_finite,
        data = x)

}#CHECK.DATA.FRAME.FINITE

# check logical flags.
check.logical = function(bool) {

  if (!is.logical(bool) || is.na(bool) || (length(bool) != 1))
    stop(sprintf("%s must be a logical value (TRUE/FALSE).",
           deparse(substitute(bool))))

}#CHECK.LOGICAL

# check a vector of weights.
check.weights = function(weights, len) {

  if (missing(weights) || is.null(weights)) {

    weights = rep(1, len)

  }#THEN
  else {

    if (!is.nonnegative.vector(weights))
      stop("missing or negative weights are not allowed.")

    if (length(weights) != len)
      stop("wrong number of weights, ", length(weights),
        " weights while ", len, " are needed.")

    weights = prop.table(weights)

  }#ELSE

  return(weights)

}#CHECK.WEIGHTS

