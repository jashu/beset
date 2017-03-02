library(beset)
context("Prediction Metrics")

test_that("lm cross-entropy = -logLik/N", {
  object <- lm(uptake ~ conc, data = CO2)
  metrics <- predict_metrics(object)
  expect_equal(metrics$mean_cross_entropy,
               -1 * as.numeric(logLik(object))/nrow(object$model))
})

test_that("gauss metrics = -logLik/N", {
  object <- glm(uptake ~ conc, data = CO2, family = "gaussian")
  metrics <- predict_metrics(object)
  expect_equal(metrics$mean_cross_entropy,
               -1 * as.numeric(logLik(object))/nrow(object$model))
})

test_that("binomial cross entropy = -logLik/N", {
  object <- glm(Treatment ~ uptake + conc, data = CO2, family = "binomial")
  metrics <- predict_metrics(object)
  expect_equal(metrics$mean_cross_entropy,
               -1 * as.numeric(logLik(object))/nrow(object$model))
})

test_that("poisson cross entropy = -logLik/N", {
  object <- glm(count ~ spray, data = InsectSprays, family = "poisson")
  metrics <- predict_metrics(object)
  expect_equal(metrics$mean_cross_entropy,
               -1 * as.numeric(logLik(object))/nrow(object$model))
})

test_that("neg binomial cross entropy = -logLik/N", {
  object <- MASS::glm.nb(count ~ spray, data = InsectSprays)
  metrics <- predict_metrics(object)
  expect_equal(metrics$mean_cross_entropy,
               -1 * as.numeric(logLik(object))/nrow(object$model))
})

test_that("zeroinfl poisson cross entropy = -logLik/N", {
  InsectSprays$count[1:10] <- 0
  object <- pscl::zeroinfl(count ~ spray, data = InsectSprays)
  metrics <- predict_metrics(object)
  expect_equal(metrics$mean_cross_entropy,
               -1 * as.numeric(logLik(object))/nrow(object$model))
})

test_that("zeroinfl negbin cross entropy = -logLik/N", {
  InsectSprays$count[1:10] <- 0
  object <- pscl::zeroinfl(count ~ spray, data = InsectSprays, dist = "negbin")
  metrics <- predict_metrics(object)
  expect_equal(metrics$mean_cross_entropy,
               -1 * as.numeric(logLik(object))/nrow(object$model))
})
