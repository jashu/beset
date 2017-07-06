library(beset)
context("Catch Linear Dependencies")

test_data <- tibble::data_frame(X1 = round(rnorm(20), 1),
                                X2 = round(rnorm(20), 1),
                                X3 = round(rnorm(20), 1))
test_data$X4 <- 0.5 * test_data$X1 -0.25 * test_data$X2 - 0.25 * test_data$X3
test_data$X5 <- factor(c(rep(1,4), rep(2,6), rep(3,10)))
test_data$X6 <- factor(c(rep(1,10), rep(2,10)))
test_data$Y <- 0.5 * test_data$X1 -0.25 * test_data$X2 - 0.25 * test_data$X3 +
  round(rnorm(20), .1)

test_that("model_frame returns same data as model.frame with .", {
  for(i in 1:7) test_data[i,i] <- NA
  mf1 <- model.frame(Y ~ ., test_data)
  mf1 <- as_tibble(mf1)
  attr(mf1, "terms") <- NULL
  attr(mf1, "na.action") <- NULL
  mf2 <- model_frame(Y ~ ., test_data)
  expect_equal(mf1, mf2)
})

test_that("model_frame returns same data as model.frame with RHS args", {
  for(i in 1:7) test_data[i,i] <- NA
  mf1 <- model.frame(Y ~ X1 + X2, test_data)
  mf1 <- as_tibble(mf1)
  attr(mf1, "terms") <- NULL
  attr(mf1, "na.action") <- NULL
  mf2 <- model_frame(Y ~ X1 + X2, test_data)
  expect_equal(mf1, mf2)
})

test_that("linear dependencies eliminated by rm_lindep", {
  mf <- model_frame(Y ~ ., test_data)
  mf_expect <- mf[c(1:4, 6)]
  mf_return <- check_lindep(Y ~ ., mf)
  expect_identical(mf_expect, mf_return)
})

test_that("rm_lindep works with non-standard variable names", {
  names(test_data) <- c("`1X`", "`X 2`", "`3`", "`X-4`", "X5", "`6 X`", "Y")
  mf <- model_frame(Y ~ ., test_data)
  mf_expect <- mf[c(1:4, 6)]
  mf_return <- check_lindep(Y ~ ., mf)
  expect_identical(mf_expect, mf_return)
  })

test_that("rm_lindep works with specified RHS for formula", {
  names(test_data) <- c("`1X`", "`X 2`", "`3`", "`X-4`", "X5", "`6 X`", "Y")
  mf <- model_frame(Y ~ `1X` + `X 2` + `3` + `X-4`, test_data)
  mf_expect <- model_frame(Y ~ `1X` + `X 2` + `3`, test_data)
  mf_return <- check_lindep(Y ~ `1X` + `X 2` + `3` + `X-4`, mf)
  expect_identical(mf_expect, mf_return)
})


