
# take care of meaningless dots arguments in plot functions.
sanitize.plot.dots = function(dots, meaningless) {

  # warn about them.
  if (any(names(dots) %in% meaningless))
    warning("arguments ", paste(meaningless, collapse = ", "),
      " will be silently ignored.")
  # nuke them from orbit.
  for (m in meaningless)
    dots[[m]] = NULL

  return(dots)

}#SANITIZE.PLOT.DOTS

# check a colour identifier (not necessarily a string/integer).
check.colour = function(col, num = 1, expand = FALSE, labels) {

  # if the colours are provided as a list (as opposed to an array) then unlist.
  if (is.list(col))
    col = unlist(col)

  # either one colour for all entities, or one colour for each entity.
  if (length(col) %!in% c(1, num)) {

    if (num == 1)
      stop(sprintf("%s must be a single colour.", deparse(substitute(col))))
    else
      stop(sprintf("%s must be a single colour or a vector of %d colours.",
        deparse(substitute(col)), num))

  }#THEN

  # this incantation returns TRUE if col2rgb() succeeds, because its return
  # value is guaranteed to be a matrix; and FALSE if col2rgb() fails and with
  # some error.
  valid = sapply(col, function(col)
              tryCatch(is.matrix(col2rgb(col)), error = function(x) { FALSE }))

  if (any(!valid))
    stop("invalid colour identifier(s) in ", deparse(substitute(col)), " :",
           paste0(" '", col[!valid], "'"))

  if (expand) {

    # expand the colour into a vector if only one colour is provided.
    if (length(col) == 1)
      col = rep(col, num)

    # add names.
    if (!missing(labels)) {

      if (is.null(names(col)))
        col = structure(col, names = labels)
      else if (!setequal(names(col), labels))
        stop("colours specified with unknown names",
          paste0(" '", col[names(col) %!in% labels], "'"))

    }#THEN

  }#THEN

  return(col)

}#CHECK.COLOUR

# check the line type identifier.
check.lty = function(lty) {

  lty.strings = c("blank", "solid", "dashed", "dotted", "dotdash", "longdash",
                  "twodash")

  if (length(lty) > 1)
    stop(sprintf("%s must be a single line type identifier.",
      deparse(substitute(lty))))
  if ((lty %!in% 0:6) && (lty %!in% lty.strings))
    stop(sprintf("%s is not a valid line type identifier.",
           deparse(substitute(lty))))

}#CHECK.LTY

# check the grid used to make probability bar/line plots easier to read.
check.quantile.grid = function(grid) {

  if (is.logical(grid)) {

    if (identical(grid, TRUE))
      grid = c(0, 0.25, 0.50, 0.75)
    else
      grid = NULL

  }#THEN
  else if (is.probability.vector(grid, zero = TRUE) && (length(grid) >= 1)) {

    if (anyDuplicated(grid))
      warning("duplicated grid points.")

  }#THEN
  else {

    stop("the grid is defined by one or more numbers between zero and one.")

  }#ELSE

  return(grid)

}#CHECK.QUANTILE.GRID
