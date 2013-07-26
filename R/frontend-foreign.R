
# read a BIF file into a bn.fit object.
read.bif = function(file, debug = FALSE) {

  # load the BIF file into memory.
  lines = readLines(file)

  read.foreign.backend(lines, format = "bif", filename = file, debug = debug)

}#READ.BIF

# dump a bn.fit object into a BIF file.
write.bif = function(file, fitted) {

  # check fitted's class.
  check.fit(fitted)
  # only discrete bayesian networks are supported.
  if (is.fitted.continuous(fitted))
    stop("only discrete Bayesian networks can be exported into BIF format.")

  # open the file for writing.
  fd = file(description = file, open = "w")

  write.foreign.backend(fd, fitted = fitted, format = "bif")

  close(fd)

  invisible(NULL)

}#WRITE.BIF

# read a DSC file into a bn.fit object.
read.dsc = function(file, debug = FALSE) {

  # load the DSC file into memory.
  lines = readLines(file)

  read.foreign.backend(lines, format = "dsc", filename = file, debug = debug)

}#READ.DSC

# dump a bn.fit object into a DSC file.
write.dsc = function(file, fitted) {

  # check fitted's class.
  check.fit(fitted)
  # only discrete bayesian networks are supported.
  if (is.fitted.continuous(fitted))
    stop("only discrete Bayesian networks can be exported into DSC format.")

  # open the file for writing.
  fd = file(description = file, open = "w")

  write.foreign.backend(fd, fitted = fitted, format = "dsc")

  close(fd)

  invisible(NULL)

}#WRITE.DSC

# read a NET file into a bn.fit object.
read.net = function(file, debug = FALSE) {

  # load the NET file into memory.
  lines = readLines(file)

  read.foreign.backend(lines, format = "net", filename = file, debug = debug)

}#READ.NET

# dump a bn.fit object into a NET file.
write.net = function(file, fitted) {

  # check fitted's class.
  check.fit(fitted)
  # only discrete bayesian networks are supported.
  if (is.fitted.continuous(fitted))
    stop("only discrete Bayesian networks can be exported into DSC format.")

  # open the file for writing.
  fd = file(description = file, open = "w")

  write.foreign.backend(fd, fitted = fitted, format = "net")

  close(fd)

  invisible(NULL)

}#WRITE.NET

