library(beset)
context("Prediction Metrics")

test_that("lm cross-entropy = -logLik/N", {
  object <- lm(uptake ~ conc, data = CO2)
  metrics <- prediction_metrics(object)
  expect_equal(metrics$mean_cross_entropy,
               -1 * as.numeric(logLik(object))/nrow(object$model))
  expect_equal(metrics$R_squared, r2d2(object))
})

test_that("gauss metrics = -logLik/N", {
  object <- glm(uptake ~ conc, data = CO2, family = "gaussian")
  metrics <- prediction_metrics(object)
  expect_equal(metrics$mean_cross_entropy,
               -1 * as.numeric(logLik(object))/nrow(object$model))
  expect_equal(metrics$R_squared, r2d2(object))
})

test_that("binomial cross entropy = -logLik/N", {
  object <- glm(Treatment ~ uptake + conc, data = CO2, family = "binomial")
  metrics <- prediction_metrics(object)
  expect_equal(metrics$mean_cross_entropy,
               -1 * as.numeric(logLik(object))/nrow(object$model))
  expect_equal(metrics$deviance_explained, r2d2(object))
})

test_that("poisson cross entropy = -logLik/N", {
  object <- glm(count ~ spray, data = InsectSprays, family = "poisson")
  metrics <- prediction_metrics(object)
  expect_equal(metrics$mean_cross_entropy,
               -1 * as.numeric(logLik(object))/nrow(object$model))
  expect_equal(metrics$deviance_explained, r2d2(object))
})

test_that("neg binomial cross entropy = -logLik/N", {
  object <- MASS::glm.nb(count ~ spray, data = InsectSprays)
  metrics <- prediction_metrics(object)
  expect_equal(metrics$mean_cross_entropy,
               -1 * as.numeric(logLik(object))/nrow(object$model))
  expect_equal(metrics$deviance_explained, r2d2(object))
})

test_that("zeroinfl poisson cross entropy = -logLik/N", {
  InsectSprays$count[1:10] <- 0
  object <- pscl::zeroinfl(count ~ spray, data = InsectSprays)
  metrics <- prediction_metrics(object)
  expect_equal(metrics$mean_cross_entropy,
               -1 * as.numeric(logLik(object))/nrow(object$model))
  expect_equal(metrics$deviance_explained, r2d2(object))
})

test_that("zeroinfl negbin cross entropy = -logLik/N", {
  InsectSprays$count[1:10] <- 0
  object <- pscl::zeroinfl(count ~ spray, data = InsectSprays, dist = "negbin")
  metrics <- prediction_metrics(object)
  expect_equal(metrics$mean_cross_entropy,
               -1 * as.numeric(logLik(object))/nrow(object$model))
  expect_equal(metrics$deviance_explained, r2d2(object))
})
