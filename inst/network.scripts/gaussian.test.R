a = data.frame(A = rnorm(5000, 1, 1),
               B = rnorm(5000, 2, 3),
               C = rep(0, 5000),
               D = rep(0, 5000),
               E = rnorm(5000, 3.5, 2),
               F = rep(0, 5000),
               G = rnorm(5000, 5, 2))
a$C = 2 * (a$A + a$B) + rnorm(5000, 2, 0.5)
a$D = 1.5 * a$B + rnorm(5000, 6, 0.33)
a$F = 2 * a$A + a$D + a$E + 1.5 * a$G + rnorm(5000, 0, 1)

