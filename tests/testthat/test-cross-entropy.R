library(beset)
context("Prediction Metrics")

test_that("lm cross-entropy = -logLik/N", {
  object <- stats::lm(Fertility ~ ., data = swiss)
  metrics <- predict_metrics(object, object$model)
  expect_equal(metrics$mce,
               -1 * as.numeric(stats::logLik(object))/nrow(object$model))
})

test_that("lm R_squared = summary()$r.squared",{
  object <- stats::lm(Fertility ~ ., data = swiss)
  metrics <- predict_metrics(object, object$model)
  expect_equal(metrics$rsq, summary(object)$r.squared)
})

test_that("gauss cross-entropy = -logLik/N", {
  object <- stats::glm(Fertility ~ ., data = swiss, family = "gaussian")
  metrics <- predict_metrics(object, object$model)
  expect_equal(metrics$mce,
               -1 * as.numeric(stats::logLik(object))/nrow(object$model))
})

test_that("gauss R_squared = 1 - deviance/null_dev",{
  object <- stats::glm(Fertility ~ ., data = swiss, family = "gaussian")
  r2 <- 1 - (object$deviance/object$null.deviance)
  metrics <- predict_metrics(object, object$model)
  expect_equal(metrics$rsq, r2)
})

test_that("binomial cross entropy = -logLik/N", {
  object <- stats::glm(Treatment ~ uptake + conc, data = CO2, family = "binomial")
  metrics <- predict_metrics(object, object$model)
  expect_equal(metrics$mce,
               -1 * as.numeric(stats::logLik(object))/nrow(object$model))
})

test_that("binomial  R_squared = 1 - deviance/null_dev)",{
  object <- stats::glm(Treatment ~ uptake + conc, data = CO2, family = "binomial")
  r2 <- 1 - (object$deviance/object$null.deviance)
  metrics <- predict_metrics(object, object$model)
  expect_equal(metrics$rsq, r2)
})

test_that("poisson cross entropy = -logLik/N", {
  object <- stats::glm(count ~ spray, data = InsectSprays, family = "poisson")
  metrics <- predict_metrics(object, object$model)
  expect_equal(metrics$mce,
               -1 * as.numeric(stats::logLik(object))/nrow(object$model))
})

test_that("poisson  R_squared = 1 - deviance/null_dev",{
  object <- stats::glm(count ~ spray, data = InsectSprays, family = "poisson")
  r2 <- 1 - (object$deviance/object$null.deviance)
  metrics <- predict_metrics(object, object$model)
  expect_equal(metrics$rsq, r2)
})

test_that("neg binomial cross entropy = -logLik/N", {
  object <- MASS::glm.nb(count ~ spray, data = InsectSprays)
  metrics <- predict_metrics(object, object$model)
  expect_equal(metrics$mce,
               -1 * as.numeric(stats::logLik(object))/nrow(object$model))
})

test_that("neg binomial  R_squared = 1 - deviance/null_dev",{
  object <- MASS::glm.nb(count ~ spray, data = InsectSprays)
  r2 <- 1 - (object$deviance/object$null.deviance)
  metrics <- predict_metrics(object, object$model)
  expect_equal(metrics$rsq, r2)
})

test_that("zeroinfl poisson cross entropy = -logLik/N", {
  object <- pscl::zeroinfl(art ~ ., data = pscl::bioChemists, dist = "poisson")
  metrics <- predict_metrics(object, object$model)
  expect_equal(metrics$mce,
               -1 * as.numeric(stats::logLik(object))/nrow(object$model))
})

test_that("zeroinfl negbin cross entropy = -logLik/N", {
  object <- pscl::zeroinfl(art ~ ., data = pscl::bioChemists, dist = "negbin")
  metrics <- predict_metrics(object, object$model)
  expect_equal(metrics$mce,
               -1 * as.numeric(stats::logLik(object))/nrow(object$model))
})
