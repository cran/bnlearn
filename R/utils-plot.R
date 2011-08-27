# print an underlined label in the plot.
underlined = function(x, y, label, col){

  text(x, y, label, col = col, font = 2)
  sw = strwidth(label)
  sh = strheight(label)
  lines(x + c(-sw/2, sw/2), rep(y - 1.5*sh/2, 2), col = col)

}#UNDERLINED

