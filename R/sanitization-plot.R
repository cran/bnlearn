
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

  # this conversion is faithful to the output of col2rgb() and prevents possible
  # "missing value where TRUE/FALSE needed" errors.
  if (any(is.na(col)))
    col[is.na(col)] = "transparent"

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
check.chart.grid = function(grid, fitted, range) {

  nodes = nodes(fitted)
  nnodes = length(nodes)
  grid.list = structure(vector(nnodes, mode = "list"), names = nodes)

  sanitize = function(grid, fitted.node, range) {

    # remove duplicated values.
    if (anyDuplicated(grid))
      warning("duplicated grid points.")
    grid = unique(grid)

    if (is(fitted.node, c("bn.fit.dnode", "bn.fit.onode"))) {

      # make sure the grid points are within the appropriate range.
      if (!is.probability.vector(grid, zero = TRUE))
        stop("grid values must be between zero and one in node ", node, ".")

    }#THEN
    else if (is(fitted.node, c("bn.fit.gnode", "bn.fit.cgnode"))) {

      # make sure there are no NAs, Infs, ets.
      if (!is.real.vector(grid))
        stop("grid values for node ", node, " must be real numbers.")

      # make sure the grid points are within the appropriate range.
      if (any(grid[grid != 0] < range[1]) || any(grid[grid != 0] > range[2]))
        warning("discarding out-of-range grid values in node ", node, ".")

      # make speciall allowances for zero, even if it is out of range.
      grid =
        union(grid[grid >= range[1] & grid <= range[2]], intersect(grid, 0))

    }#THEN

  }#SANITIZE

  if (is.logical(grid)) {

    for (node in nodes) {

      if (is(fitted[[node]], c("bn.fit.dnode", "bn.fit.onode"))) {

        if (identical(grid, TRUE))
          grid.list[[node]] = c(0, 0.25, 0.50, 0.75)
        else
          grid.list[node] = list(NULL)

      }#THEN
      else if (is(fitted[[node]], c("bn.fit.gnode", "bn.fit.cgnode"))) {

        if (identical(grid, TRUE)) {

          if (all(range[[node]] <= 0) || all(range[[node]] >= 0))
            grid.list[[node]] = c(0, range[[node]], mean(range[[node]]))
          else
            grid.list[[node]] = c(0, range[[node]] / 2, range[[node]])

        }#THEN
        else {

          grid.list[node] = list(NULL)

        }#ELSE

      }#THEN

    }#FOR

  }#THEN
  else if (is.real.vector(grid) && (length(grid) >= 1)) {

    for (node in nodes)
      grid.list[[node]] = sanitize(grid, fitted[[node]], range[[node]])

  }#THEN
  else if (is.list(grid)) {

    if (length(grid) != nnodes)
      stop("the list of grid points has ", length(grid),
           " elements, but the network has ", nnodes, " nodes.")
    if (is.null(names(grid)) || !setequal(names(grid), nodes))
      stop("the list of grid points must be a named list with one element per node.")

    for (node in nodes)
      grid.list[[node]] = sanitize(grid[[node]], fitted[[node]], range[[node]])

  }#THEN
  else {

    stop("the grid is defined by one or more numbers in the range of the parameter values.")

  }#ELSE

  return(grid.list)

}#CHECK.CHART.GRID
