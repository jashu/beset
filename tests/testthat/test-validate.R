library(beset)
context("validate")

test_that("lm validate yields plausible rsq", {
  object <- stats::lm(Fertility ~ ., data = swiss)
  results <- validate(object)
  expect_true(between(results$stats$rsq$mean, 0.58, 0.62))
})

test_that("glm validate yields plausible rsq", {
  object <- stats::glm(Treatment ~ uptake + conc, data = CO2,
                       family = "binomial")
  results <- validate(object)
  expect_true(between(results$stats$rsq$mean, 0.04, 0.07))
})

test_that("neg binomial validate yields plausible rsq",{
  object <- MASS::glm.nb(count ~ spray, data = InsectSprays)
  results <- validate(object)
  expect_true(between(results$stats$rsq$mean, 0.7, 0.73))
})

object <- pscl::zeroinfl(art ~ ., data = pscl::bioChemists, dist = "poisson")
test_that(
  "internal zi.fit returns same coefs as zeroinfl for Poisson family", {
    zif_par <- beset:::set_zi_par(object)
    zif_fit <- do.call(beset:::zi.fit, zif_par)
    expect_identical(object$coefficients, zif_fit$coefficients)
  })
test_that("zeroinfl poisson validate yields plausible rsq", {
  results <- validate(object)
  expect_true(between(results$stats$rsq$mean, 0.12, 0.14))
})

object <- pscl::zeroinfl(art ~ ., data = pscl::bioChemists, dist = "negbin")
test_that("internal zi.fit returns same coefs as zeroinfl for negbin family", {
  zif_par <- beset:::set_zi_par(object)
  zif_fit <- do.call(beset:::zi.fit, zif_par)
  expect_identical(object$coefficients, zif_fit$coefficients)
})

test_that("zeroinfl negbin validate yields plausible rsq", {
  results <- validate(object)
  expect_true(between(results$stats$rsq$mean, 0.09, 0.1))
})

