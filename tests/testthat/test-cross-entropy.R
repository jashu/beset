library(beset)
context("Cross entropy")

test_that("gauss cross entropy = -logLik/N", {
  object <- glm(uptake ~ conc, data = CO2, family = "gaussian")
  expect_equal(cross_entropy(object, object$model),
               -1 * as.numeric(logLik(object))/nrow(object$model))
})

test_that("binomial cross entropy = -logLik/N", {
  object <- glm(Treatment ~ uptake + conc, data = CO2, family = "binomial")
  expect_equal(cross_entropy(object, object$model),
               -1 * as.numeric(logLik(object))/nrow(object$model))
})

test_that("poisson cross entropy = -logLik/N", {
  object <- glm(count ~ spray, data = InsectSprays, family = "poisson")
  expect_equal(cross_entropy(object, object$model),
               -1 * as.numeric(logLik(object))/nrow(object$model))
})

test_that("neg binomial cross entropy = -logLik/N", {
  object <- MASS::glm.nb(count ~ spray, data = InsectSprays)
  expect_equal(cross_entropy(object, object$model),
               -1 * as.numeric(logLik(object))/nrow(object$model))
})

test_that("zeroinfl poisson cross entropy = -logLik/N", {
  InsectSprays$count[1:10] <- 0
  object <- pscl::zeroinfl(count ~ spray, data = InsectSprays)
  expect_equal(cross_entropy(object, object$model),
               -1 * as.numeric(logLik(object))/nrow(object$model))
})

test_that("zeroinfl negbin cross entropy = -logLik/N", {
  InsectSprays$count[1:10] <- 0
  object <- pscl::zeroinfl(count ~ spray, data = InsectSprays, dist = "negbin")
  expect_equal(cross_entropy(object, object$model),
               -1 * as.numeric(logLik(object))/nrow(object$model))
})
