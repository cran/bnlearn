a = sample(c("a", "b", "c"), 5000, prob = rep(1/3, 3), replace = TRUE)
c = sample(c("a", "b", "c"), 5000, prob = c(0.75, 0.2, 0.05),
      replace = TRUE)
f = sample(c("a", "b"), 5000, prob = rep(1/2, 2), replace = TRUE)

b = a
b[b == "a"] = sample(c("a", "b", "c"), length(which(b == "a")),
                prob = c(0.8, 0.1, 0.1), replace = TRUE)
b[b == "b"] = sample(c("a", "b", "c"), length(which(b == "b")),
                prob = c(0.4, 0.2, 0.4), replace = TRUE)
b[b == "c"] = sample(c("a", "b", "c"), length(which(b == "c")),
                prob = c(0.1, 0.1, 0.8), replace = TRUE)

d = apply(cbind(a,c), 1, paste, collapse= ":")
d[d == "a:a"] = sample(c("a", "b", "c"), length(which(d == "a:a")),
                  prob = c(0.8, 0.1, 0.1), replace = TRUE)
d[d == "a:b"] = sample(c("a", "b", "c"), length(which(d == "a:b")),
                  prob = c(0.2, 0.1, 0.7), replace = TRUE)
d[d == "a:c"] = sample(c("a", "b", "c"), length(which(d == "a:c")),
                  prob = c(0.4, 0.2, 0.4), replace = TRUE)
d[d == "b:a"] = sample(c("a", "b", "c"), length(which(d == "b:a")),
                  prob = c(0.1, 0.8, 0.1), replace = TRUE)
d[d == "b:b"] = sample(c("a", "b", "c"), length(which(d == "b:b")),
                  prob = c(0.9, 0.05, 0.05), replace = TRUE)
d[d == "b:c"] = sample(c("a", "b", "c"), length(which(d == "b:c")),
                  prob = c(0.3, 0.4, 0.3), replace = TRUE)
d[d == "c:a"] = sample(c("a", "b", "c"), length(which(d == "c:a")),
                  prob = c(0.1, 0.1, 0.8), replace = TRUE)
d[d == "c:b"] = sample(c("a", "b", "c"), length(which(d == "c:b")),
                  prob = c(0.25, 0.5, 0.25), replace = TRUE)
d[d == "c:c"] = sample(c("a", "b", "c"), length(which(d == "c:c")),
                  prob = c(0.15, 0.45, 0.4), replace = TRUE)

e = apply(cbind(b,f), 1, paste, collapse= ":")
e[e == "a:a"] = sample(c("a", "b", "c"), length(which(e == "a:a")),
                  prob = c(0.8, 0.1, 0.1), replace = TRUE)
e[e == "a:b"] = sample(c("a", "b", "c"), length(which(e == "a:b")),
                  prob = c(0.4, 0.5, 0.1), replace = TRUE)
e[e == "b:a"] = sample(c("a", "b", "c"), length(which(e == "b:a")),
                  prob = c(0.2, 0.2, 0.6), replace = TRUE)
e[e == "b:b"] = sample(c("a", "b", "c"), length(which(e == "b:b")),
                  prob = c(0.3, 0.4, 0.3), replace = TRUE)
e[e == "c:a"] = sample(c("a", "b", "c"), length(which(e == "c:a")),
                  prob = c(0.1, 0.1, 0.8), replace = TRUE)
e[e == "c:b"] = sample(c("a", "b", "c"), length(which(e == "c:b")),
                  prob = c(0.25, 0.5, 0.25), replace = TRUE)

data.frame(A = a, B = b, C = c, D = d, E = e, F = f)

