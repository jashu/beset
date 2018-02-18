library(beset)
context("Prediction Metrics")

test_that("lm cross-entropy = -logLik/N", {
  object <- stats::lm(Fertility ~ ., data = swiss)
  metrics <- predict_metrics(object, object$model)
  expect_equal(metrics$mean_cross_entropy,
               -1 * as.numeric(stats::logLik(object))/nrow(object$model))
})

test_that("lm R_squared = summary()$r.squared",{
  object <- stats::lm(Fertility ~ ., data = swiss)
  metrics <- predict_metrics(object, object$model)
  expect_equal(metrics$R_squared, summary(object)$r.squared)
})

test_that("lm deviance = deviance()",{
  object <- stats::lm(Fertility ~ ., data = swiss)
  metrics <- predict_metrics(object, object$model)
  expect_equal(metrics$deviance, stats::deviance(object))
})

test_that("gauss cross-entropy = -logLik/N", {
  object <- stats::glm(Fertility ~ ., data = swiss, family = "gaussian")
  metrics <- predict_metrics(object, object$model)
  expect_equal(metrics$mean_cross_entropy,
               -1 * as.numeric(stats::logLik(object))/nrow(object$model))
})

test_that("gauss deviance = deviance()",{
  object <- stats::glm(Fertility ~ ., data = swiss, family = "gaussian")
  metrics <- predict_metrics(object, object$model)
  expect_equal(metrics$deviance, stats::deviance(object))
})

test_that("binomial cross entropy = -logLik/N", {
  object <- stats::glm(Treatment ~ uptake + conc, data = CO2, family = "binomial")
  metrics <- predict_metrics(object, object$model)
  expect_equal(metrics$mean_cross_entropy,
               -1 * as.numeric(stats::logLik(object))/nrow(object$model))
})

test_that("binomial deviance = deviance()",{
  object <- stats::glm(Treatment ~ uptake + conc, data = CO2, family = "binomial")
  metrics <- predict_metrics(object, object$model)
  expect_equal(metrics$deviance, stats::deviance(object))
})

test_that("poisson cross entropy = -logLik/N", {
  object <- stats::glm(count ~ spray, data = InsectSprays, family = "poisson")
  metrics <- predict_metrics(object, object$model)
  expect_equal(metrics$mean_cross_entropy,
               -1 * as.numeric(stats::logLik(object))/nrow(object$model))
})

test_that("poisson deviance = deviance()",{
  object <- stats::glm(count ~ spray, data = InsectSprays, family = "poisson")
  metrics <- predict_metrics(object, object$model)
  expect_equal(metrics$deviance, stats::deviance(object))
})

test_that("neg binomial cross entropy = -logLik/N", {
  object <- MASS::glm.nb(count ~ spray, data = InsectSprays)
  metrics <- predict_metrics(object, object$model)
  expect_equal(metrics$mean_cross_entropy,
               -1 * as.numeric(stats::logLik(object))/nrow(object$model))
})

test_that("neg binomial deviance = deviance()",{
  object <- MASS::glm.nb(count ~ spray, data = InsectSprays)
  metrics <- predict_metrics(object, object$model)
  expect_equal(metrics$deviance, stats::deviance(object))
})

test_that("zeroinfl poisson cross entropy = -logLik/N", {
  object <- pscl::zeroinfl(art ~ ., data = pscl::bioChemists, dist = "poisson")
  metrics <- predict_metrics(object, object$model)
  expect_equal(metrics$mean_cross_entropy,
               -1 * as.numeric(stats::logLik(object))/nrow(object$model))
})

test_that("zeroinfl negbin cross entropy = -logLik/N", {
  object <- pscl::zeroinfl(art ~ ., data = pscl::bioChemists,
                           dist = "negbin")
  metrics <- predict_metrics(object, object$model)
  expect_equal(metrics$mean_cross_entropy,
               -1 * as.numeric(stats::logLik(object))/nrow(object$model))
})

test_that("zeroinfl poisson cross entropy = -logLik/N", {
  data(BinomialExample, package = "glmnet")
  object <- glmnet::glmnet(x, y, family = "binomial")
  metrics <- predict_metrics(object, object$model)
  expect_equal(metrics$mean_cross_entropy,
               -1 * as.numeric(stats::logLik(object))/nrow(object$model))
})

test_that("zeroinfl negbin cross entropy = -logLik/N", {
  object <- pscl::zeroinfl(art ~ ., data = pscl::bioChemists,
                           dist = "negbin")
  metrics <- predict_metrics(object, object$model)
  expect_equal(metrics$mean_cross_entropy,
               -1 * as.numeric(stats::logLik(object))/nrow(object$model))
})
