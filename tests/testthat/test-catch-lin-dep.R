library(beset)
context("Catch Linear Dependencies")

test_that("linear dependencies eliminated by rm_lindep", {
  test_data <- tibble::data_frame(X1 = round(rnorm(20), 1),
                                  X2 = round(rnorm(20), 1),
                                  X3 = round(rnorm(20), 1))
  test_data$X4 <- 0.5 * test_data$X1 -0.25 * test_data$X2 - 0.25 * test_data$X3
  test_data$X5 <- factor(c(rep(1,4), rep(2,6), rep(3,10)))
  test_data$X6 <- factor(c(rep(1,10), rep(2,10)))
  test_data$Y <- 0.5 * test_data$X1 -0.25 * test_data$X2 - 0.25 * test_data$X3 +
    round(rnorm(20), .1)
  X <- model.frame(Y ~ ., test_data)
  lin_dep <- which(colnames(X) %in% c("X4", "X6"))
  X_expect <- X[, -lin_dep]
  X_return <- beset:::check_lindep(X)
  expect_identical(X_expect, X_return)
})

