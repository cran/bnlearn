library(bnlearn)

underflow = data.frame(
  A = factor("a", levels=c("a", "b")),
  B = sample(c(rep("a", 12), rep("b", 188)), 200))
table(underflow)
format(ci.test(x = "A", y = "B", data = underflow, test = "mi")$statistic, digits = 20)
format(ci.test(x = "B", y = "A", data = underflow, test = "mi")$statistic, digits = 20)
overflow = as.data.frame(matrix(c(rep(c("A", "A"), 32114), rep(c("A", "B"), 40677), rep(c("B", "A"), 77472), rep(c("B", "B"), 126716)), ncol = 2, byrow = TRUE))
overflow[, 1] = as.factor(overflow[, 1])
overflow[, 2] = as.factor(overflow[, 2])
names(overflow) = LETTERS[1:2]
table(overflow)
format(ci.test(x = "A", y = "B", data = overflow, test = "mi")$statistic, digits = 20)
format(ci.test(x = "B", y = "A", data = overflow, test = "mi")$statistic, digits = 20)

