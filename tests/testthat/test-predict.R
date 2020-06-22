# library(beset)
# context("predict")
#
# test_that("predict method for beset_lm works", {
#   object <- beset_lm(Fertility ~ ., data = swiss)
#   yhat <- predict(object)
#   expect_length(yhat, 47)
# })
#
# test_that("predict method for beset_lm works with nested cv", {
#   object <- beset_lm(Fertility ~ ., data = swiss, nest_cv = TRUE)
#   yhat <- predict(object)
#   expect_length(yhat, 47)
# })
#
# test_that("predict method for beset_elnet works", {
#   object <- beset_elnet(Fertility ~ ., data = swiss)
#   yhat <- predict(object)
#   expect_length(yhat, 47)
# })
#
# test_that("predict method for beset_elnet works with nested cv", {
#   object <- beset_elnet(Fertility ~ ., data = swiss, nest_cv = TRUE)
#   yhat <- predict(object)
#   expect_length(yhat, 47)
# })
#
# test_that("predict method for glm binomial family yields plausible numbers", {
#   object <- beset_glm(Treatment ~ uptake + conc, data = CO2, family = "binomial")
#   yhat <- predict(object)
#   expect_true(all(between(yhat, 0, 1)))
# })
#
# test_that("predict method for elnet binomial family yields plausible numbers", {
#     object <- beset_elnet(Treatment ~ ., data = CO2, family = "binomial")
#     yhat <- predict(object)
#     expect_true(all(between(yhat, 0, 1)))
# })
