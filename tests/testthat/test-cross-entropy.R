library(beset)
context("Prediction Metrics")

test_that("lm cross-entropy = -logLik/N and R_squared = summary()$r.squared", {
  object <- lm(Fertility ~ ., data = swiss)
  metrics <- predict_metrics(object)
  expect_equal(metrics$mean_cross_entropy,
               -1 * as.numeric(logLik(object))/nrow(object$model))
  expect_equal(metrics$R_squared, summary(object)$r.squared)
})

test_that("gauss metrics = -logLik/N", {
  object <- glm(Fertility ~ ., data = swiss, family = "gaussian")
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
  object <- pscl::zeroinfl(art ~ ., data = pscl::bioChemists,
                           dist = "negbin")
  metrics <- predict_metrics(object)
  expect_equal(metrics$mean_cross_entropy,
               -1 * as.numeric(logLik(object))/nrow(object$model))
})
