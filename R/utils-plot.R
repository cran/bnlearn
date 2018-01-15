# print an underlined label in the plot.
underlined = function(x, y, label, col){

  text(x, y, label, col = col, font = 2)
  sw = strwidth(label)
  sh = strheight(label)
  lines(x + c(-sw/2, sw/2), rep(y - 1.5*sh/2, 2), col = col)

}#UNDERLINED

# compute the largest possible expansion factor for which a string fits the box.
largest.cex = function(node, height, width, hfrac = 0.7, wfrac = 0.9) {

  # fit vertically.
  guess  = hfrac * height / strheight(node, cex = 1)
  best = optimize(f = function(x) abs(strheight(node, cex = x) - hfrac * height),
           interval = guess * c(0.5, 2), tol = 0.025)$minimum

  # fit horizonally.
  best = best * min(wfrac * width / strwidth(node, cex = best), 1)

  return(best)

}#LARGEST.CEX

# created a lighter tint of a given colour.
lighter.colour = function(col, offset = 0.25) {

  rgb = col2rgb(col)[, 1]
  do.call("rgb", as.list((rgb + (255 - rgb) * offset) / 255))

}#LIGHTER.COLOUR
