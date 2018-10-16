# advanced cat which correctly handles ini-like lines and short line widths.
wcat = function(header, value) {

  # get the number of available columns.
  columns = options("width")

  # blatantly ignore any line width shorter than 45, trying to support
  # that case is a losing proposition.
  if ((columns >= nchar(header) + nchar(value) + 1) || (columns < 45)) {

    # if there are enough columns print the string as is.
    cat(paste(header, value, sep = " "), "\n")

  }#THEN
  else {

    # if there are not enough columns print the header on one row
    # (left-aligned) and the value on the following line (right-aligned).
    cat(header, "\n", sprintf(paste("%", columns, "s", sep = ""), value), "\n")

  }#ELSE

}#WCAT

# advanced cat handles model strings and short line widths.
fcat = function(modelstr) {

  # measure the number of available columns.
  columns = options("width")

  if ((columns >= nchar(modelstr)) || columns < 45) {

    cat("  ", modelstr, "\n")

  }#THEN
  else {

    cat(paste(strsplit(modelstr, "\\]")[[1]], "]", sep = ""),
      fill = TRUE, sep = "", labels = "  ")

  }#ELSE

}#FCAT
