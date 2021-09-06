
# compute the sample size / CPT cells ratio.
obs.per.cell = function(x, y, z = NULL, data) {

  opc = 0
  ndata = nrow(data)
  nlx = nlevels(data[, x])
  nly = nlevels(data[, y])

  # return +Inf for continuous data to bypass countermeasures
  # thought for sparce discrete data.
  if (data.type(data) == "continuous")
    return(Inf)

  if (is.null(z) || (length(z) == 0)) {

    opc = ndata / (nlx * nly)

  }#THEN
  else if (is.character(z)) {

    if (length(z) == 1)
      opc = ndata / (nlx * nly * nlevels(data[, z]))
    else if (length(z) > 1)
      opc = ndata / (nlx * nly *
        prod(sapply(z, function(col) { nlevels(data[, col]) } )))

  }#THEN
  else if (is.factor(z)) {

    opc = ndata / (nlx * nly * nlevels(z))

  }#ELSE

  return(opc)

}#OBS.PER.CELL

