
# check labels for various arguments.
check.label = function(arg, choices, labels, argname, see) {

  if (!is.string(arg))
    stop("the ", argname, " must be a single character string.")

  if (arg %in% choices)
    return(invisible(NULL))

  # concatenate valid values, optinally with labels.
  if (missing(labels)) {

    choices = paste(paste('"', choices, '"', sep = ""), collapse = ", ")

  }#THEN
  else {

    labels = paste("(", labels[choices], ")", sep = "")
    choices = paste('"', choices, '"', sep = "")
    nl = length(labels)
    choices = paste(choices, labels, collapse = ", ")

  }#THEN

  # mention the most relevant manual page.
  if (missing(see))
    see = character(0)
  else
    see = paste(" See ?", see, " for details.", sep = "")

  # build the error message.
  errmsg = paste("valid ", argname, "(s) are ", choices, ".", see, sep = "")

  # make sure the error message is not truncated if possible at all.
  errlen = unlist(options("warning.length"), use.names = FALSE)
  options("warning.length" = max(1000, min(8170, nchar(errmsg) + 20)))
  on.exit(options("warning.length" = errlen))

  stop(errmsg)

}#CHECK.LABEL

# warn about unused arguments.
check.unused.args = function(dots, used.args) {

  if (is(dots, "list"))
    args = names(dots)
  else
    args = dots

  unused = args %!in% used.args

  if (any(unused))
    warning("unused argument(s):", paste0(" '", args[unused], "'"), ".")

}#CHECK.UNUSED.ARGS

# check whether a package is loaded.
check.and.load.package = function(pkg) {

  # silence all warnings while looking for suggested packages.
  warning.level  = as.numeric(options("warn"))
  options("warn" = -1)
  on.exit(options("warn" = warning.level))

  if (!requireNamespace(pkg))
    stop("this function requires the ", pkg, " package.")

}#CHECK.AND.LOAD.PACKAGE

# reverse lookup for optional arguments.
has.argument = function(label, arg, lookup) {

  arg %in% lookup[[label]]

}#HAS.ARGUMENT
