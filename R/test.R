
conditional.test = function(x, y, sx, data, test) {

  assign("test.counter", get("test.counter", envir = .GlobalEnv) + 1, 
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

      pchisq(mi.test(table(data[,x], data[,y]), gsquare = TRUE), 
        (nlevels(data[,x]) - 1) * (nlevels(data[,y]) - 1), lower.tail = FALSE)

    }#THEN
    # Mutual Infomation a la FastIAMB (chi-square asymptotic distribution)
    else if (test == "fmi") {

      if (obs.per.cell(x, y, data = data) >= 5) {

        pchisq(mi.test(table(data[,x], data[,y]), gsquare = TRUE), 
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

    config = factor(apply(as.data.frame(data[,sx]), 1, paste, 
      sep = "", collapse = ":"))

    # Cochran-Mantel-Haenszel (chi-square asymptotic distribution)
    if (test == "mh") {

      mantelhaen.test(data[,x], data[,y], config, exact=TRUE)$p.value

    }#THEN
    # Conditional Mutual Infomation (chi-square asymptotic distribution)
    else if (test == "mi") {

      pchisq(cmi.test(data[,x], data[,y], config, gsquare = TRUE), 
        (nlevels(data[,x]) - 1) * (nlevels(data[,y]) - 1) * nlevels(config), lower.tail = FALSE)

    }#THEN
    # Conditional Mutual Infomation a la FastIAMB (chi-square asymptotic distribution)
    else if (test == "fmi") {

      if (obs.per.cell(x, y, config, data = data) >= 5) {

        pchisq(cmi.test(data[,x], data[,y], config, gsquare = TRUE), 
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

mi.test = function(table, gsquare = TRUE) {

  result = 0

  s = .C("mi",
    n = table,
    nc = as.integer(dim(table)[1]),
    nr = as.integer(dim(table)[2]),
    result = as.double(result),
    PACKAGE = "bnlearn")

  if (gsquare) 
    2 * sum(table) * s$result
  else 
    s$result

}#MI.TEST

cmi.test = function(x, y, z, gsquare = TRUE) {

  result = 0
  s = .C("cmi",
      x = x,
      y = y,
      z = z,
      lx = as.integer(nlevels(x)),
      ly = as.integer(nlevels(y)),
      lz = as.integer(nlevels(z)),
      length = as.integer(length(x)),
      result = as.double(result),
      PACKAGE = "bnlearn")

  if (gsquare) 
    2 * length(x) * s$result
  else 
    s$result

}#CMI.TEST

# Partial Canonical (Linear) Correlation.
# Internal copy of the pcor() function from package ggm.
pcor = function (u, S) {

  k = solve(cov(S[, u]))
  -k[1, 2]/sqrt(k[1, 1] * k[2, 2])

}#PCOR

