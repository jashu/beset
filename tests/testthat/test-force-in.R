library(beset)
# context("forced_cols")

mf <- model.frame(~ count + spray, data = InsectSprays)
x <- model.matrix(~ count + spray, data = InsectSprays)

test_that("index of 1 is returned when force_in is NULL", {
  expect_equivalent(forced_cols(mf, x), 1)
})

test_that("indicies of numerics match", {
 expect_equivalent(forced_cols(mf, x, "count"), c(1,2))
})

test_that("indicies of factors match",{
  expect_equivalent(forced_cols(mf, x, "spray"), c(1, 3:7))
})

test_that("indicies of numerics + factors match",{
  expect_equivalent(forced_cols(mf, x, c("count", "spray")), 1:7)
})

x <-model.matrix(~ 0 + count + spray, data = InsectSprays)

test_that("indicies match when there is no intercept",{
  expect_equivalent(forced_cols(mf, x, c("count", "spray")), 1:7)
})

test_that("returns empty vector when no intercept and force_in is NULL",{
  expect_equivalent(forced_cols(mf, x), integer(0))
})
