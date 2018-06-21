library(beset)
context("glm_nb")

x <- model.matrix(~ spray, data = InsectSprays)
y <- InsectSprays$count
object_beset <- glm_nb(x, y)
object_mass <- MASS::glm.nb(count ~ spray, data = InsectSprays)

test_that("coefficients of glm.nb and glm_nb match", {
 expect_equivalent(object_mass$coefficients, object_beset$coefficients)
})

test_that("residuals of glm.nb and glm_nb match",{
  expect_equivalent(object_mass$residuals, object_beset$residuals)
})

test_that("fitted.values of glm.nb and glm_nb match",{
  expect_equivalent(object_mass$fitted.values, object_beset$fitted.values)
})

test_that("aic of glm.nb and glm_nb match",{
  expect_equivalent(object_mass$aic, object_beset$aic)
})

test_that("null.deviance of glm.nb and glm_nb match",{
  expect_equivalent(object_mass$null.deviance, object_beset$null.deviance)
})

test_that("deviance of glm.nb and glm_nb match",{
  expect_equivalent(object_mass$deviance, object_beset$deviance)
})
