
# check the string representation (aka the formula) of a network.
check.modelstring = function(string) {

  # check the type.
  if (!is.string(string))
    stop("string must be a character string.")

  # check the syntax (separate regexps for root nodes and non-root ndoes).
  correct.format = paste("^(",
    "\\[[^\\[\\]\\|:]+?\\]",
    "|",
    "\\[[^\\[\\]\\|:]+?\\|[^\\[\\]\\|:]+?([:]{0,1}[^\\[\\]\\|:])*?\\]",
  ")+$", sep = "")

  if (!grepl(correct.format, string, perl = TRUE))
    stop("malformed model string format (see ?modelstring).")

}#CHECK.MODELSTRING

