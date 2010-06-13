BOOL = c("yes", "no")

a = sample(BOOL, 5000, prob = c(0.01, 0.99), replace = TRUE)
s = sample(BOOL, 5000, prob = c(0.50, 0.50), replace = TRUE)

t = a
t[t == "yes"] = sample(BOOL, length(which(t == "yes")), prob = c(0.05, 0.95), replace = TRUE)
t[t == "no"] = sample(BOOL, length(which(t == "no")), prob = c(0.01, 0.99), replace = TRUE)

l = s
l[l == "yes"] = sample(BOOL, length(which(l == "yes")), prob = c(0.10, 0.90), replace = TRUE)
l[l == "no"] = sample(BOOL, length(which(l == "no")), prob = c(0.01, 0.99), replace = TRUE)

b = s
b[b == "yes"] = sample(BOOL, length(which(b == "yes")), prob = c(0.60, 0.40), replace = TRUE)
b[b == "no"] = sample(BOOL, length(which(b == "no")), prob = c(0.30, 0.70), replace = TRUE)

e = apply(cbind(l,t), 1, paste, collapse= ":")
e[e == "yes:yes"] = "yes"
e[e == "yes:no"] = "yes"
e[e == "no:yes"] = "yes"
e[e == "no:no"] = "no"

x = e
x[x == "yes"] = sample(BOOL, length(which(x == "yes")), prob = c(0.98, 0.02), replace = TRUE)
x[x == "no"] = sample(BOOL, length(which(x == "no")), prob = c(0.05, 0.95), replace = TRUE)

d = apply(cbind(e,b), 1, paste, collapse= ":")
d[d == "yes:yes"] = sample(BOOL, length(which(d == "yes:yes")), prob = c(0.90, 0.10), replace = TRUE)
d[d == "yes:no"] = sample(BOOL, length(which(d == "yes:no")), prob = c(0.70, 0.30), replace = TRUE)
d[d == "no:yes"] = sample(BOOL, length(which(d == "no:yes")), prob = c(0.80, 0.20), replace = TRUE)
d[d == "no:no"] = sample(BOOL, length(which(d == "no:no")), prob = c(0.10, 0.90), replace = TRUE)

asia = data.frame(
  A = factor(a, levels = BOOL),
  S = factor(s, levels = BOOL),
  T = factor(t, levels = BOOL),
  L = factor(l, levels = BOOL),
  B = factor(b, levels = BOOL),
  E = factor(e, levels = BOOL),
  X = factor(x, levels = BOOL),
  D = factor(d, levels = BOOL)
)
