
conditional.test = function(x, y, sx, data, test) {

  # update the test counter.
  assign(".test.counter", get(".test.counter", envir = .GlobalEnv) + 1,
    envir = .GlobalEnv)

  sx = sx[sx != ""]
  ndata = nrow(data)

  if (length(sx) == 0) {

    # Cochran-Mantel-Haenszel with dummy stratification (chi-square asymptotic distribution)
    if (test == "mh") {

      mantelhaen.test(data[,x], data[,y], factor(rep(1, ndata)), exact=TRUE)$p.value

    }#THEN
    # Mutual Infomation (chi-square asymptotic distribution)
    else if (test == "mi") {

      pchisq(mi.test(data[,x], data[,y], ndata, gsquare = TRUE),
        (nlevels(data[,x]) - 1) * (nlevels(data[,y]) - 1), lower.tail = FALSE)

    }#THEN
    # Akaike Information Criterion-like test (binary, no p-value!)
    else if (test == "aict") {

      as.integer(mi.test(data[,x], data[,y], ndata, gsquare = FALSE) <
        (nlevels(data[,x]) - 1) * (nlevels(data[,y]) - 1) / ndata)

    }#THEN
    # Mutual Infomation a la FastIAMB (chi-square asymptotic distribution)
    else if (test == "fmi") {

      if (obs.per.cell(x, y, data = data) >= 5) {

        pchisq(mi.test(data[,x], data[,y], ndata, gsquare = TRUE),
          (nlevels(data[,x]) - 1) * (nlevels(data[,y]) - 1), lower.tail = FALSE)

      }#THEN
      else {

        return(1);

      }#ELSE

    }#THEN
    # Canonical (Linear) Correlation (Student's t distribution)
    else if (test == "cor") {

     rxy = cor(data[,x], data[,y])
     pt(abs((rxy * sqrt(ndata - 2) / sqrt(1 - rxy^2))), ndata - 2, lower.tail = FALSE) * 2

    }#THEN
    # Fisher's Z (asymptotic normal distribution)
    else if (test == "zf") {

      rxy = cor(data[,x], data[,y])
      pnorm(abs(log((1 + rxy)/(1 - rxy))/2 * sqrt(ndata -3)), lower.tail = FALSE) * 2

    }#THEN

  }#THEN
  else {

    # build the contingency table for discrete data only.
    if (test %in% available.discrete.tests) {

      # if there is only one parent, get it easy.
      if (length(sx) == 1)
        config = data[, sx]
      else
        config = configurations(data[, sx])

    }#THEN

    # Cochran-Mantel-Haenszel (chi-square asymptotic distribution)
    if (test == "mh") {

      mantelhaen.test(data[,x], data[,y], config, exact=TRUE)$p.value

    }#THEN
    # Conditional Mutual Infomation (chi-square asymptotic distribution)
    else if (test == "mi") {

      pchisq(cmi.test(data[,x], data[,y], config, ndata, gsquare = TRUE),
        (nlevels(data[,x]) - 1) * (nlevels(data[,y]) - 1) * nlevels(config), lower.tail = FALSE)

    }#THEN
    # Conditional Akaike Information Criterion-like test (binary, no p-value!)
    else if (test == "aict") {

      as.integer(cmi.test(data[,x], data[,y], config, ndata, gsquare = FALSE) <
        (nlevels(data[,x]) - 1) * (nlevels(data[,y]) - 1) * nlevels(config) / ndata)

    }#THEN
    # Conditional Mutual Infomation a la FastIAMB (chi-square asymptotic distribution)
    else if (test == "fmi") {

      if (obs.per.cell(x, y, config, data = data) >= 5) {

        pchisq(cmi.test(data[,x], data[,y], config, ndata, gsquare = TRUE),
          (nlevels(data[,x]) - 1) * (nlevels(data[,y]) - 1) * nlevels(config), lower.tail = FALSE)

      }#THEN
      else {

        return(1);

      }#ELSE

    }#THEN
    # Canonical Partial Correlation (Student's t distribution)
    else if (test == "cor") {

      rxy.z = pcor(c(x, y, sx), data)
      df = ndata - 2 - length(sx)
      pt(abs(rxy.z * sqrt(df) / sqrt(1 - rxy.z^2)), df, lower.tail = FALSE) * 2

    }#THEN
    # Fisher's Z (asymptotic normal distribution)
    else if (test == "zf") {

      rxy.z = pcor(c(x, y, sx), data)
      df = ndata - 3 - length(sx)
      pnorm(abs(log((1 + rxy.z)/(1 - rxy.z))/2 * sqrt(df)), lower.tail = FALSE) * 2

    }#THEN

  }#ELSE

}#CONDITIONAL.TEST

mi.test = function(x, y, ndata, gsquare = TRUE) {

  s = .Call("mi",
      x = x,
      y = y,
      lx = nlevels(x),
      ly = nlevels(y),
      length = ndata,
      PACKAGE = "bnlearn")

  ifelse(gsquare, 2 * ndata * s, s)

}#MI.TEST

cmi.test = function(x, y, z, ndata, gsquare = TRUE) {

  s = .Call("cmi",
      x = x,
      y = y,
      z = z,
      lx = nlevels(x),
      ly = nlevels(y),
      lz = nlevels(z),
      length = ndata,
      PACKAGE = "bnlearn")

  ifelse(gsquare, 2 * ndata * s, s)

}#CMI.TEST

# Partial Canonical (Linear) Correlation.
# Internal copy of the pcor() function from package ggm.
# Copyright Giovanni M. Marchetti, 2006.
# Released under "GPL version 2 or newer".
pcor = function (u, S) {

  k = solve(cov(S[, u]))
  -k[1, 2]/sqrt(k[1, 1] * k[2, 2])

}#PCOR

