#' Prediction Metrics
#'
#' Calculate mean squared error, mean cross entropy, and R-squared (as fraction
#' of deviance explained) to measure how well a generalized linear model
#' predicts test data.
#'
#' \code{predict_metrics} uses a generalized linear model (GLM) previously
#' fit to a training set to predict responses for a test set. It then computes
#' the residual sum of squares and the log-likelihood of the predicted
#' responses, as well as the log-likelihoods for the corresponding saturated
#' model (a model with one free parameter per observation) and null model (a
#' model with only an intercept) fit to the true responses. These quantities
#' are used to derive the metrics described below.
#'
#' @section Mean squared error (MSE) and mean cross entropy (MCE):
#' MSE is the average of the squared difference between each predicted response
#' \eqn{ŷ} and the actual observed response \eqn{y}. MCE is an analagous
#' information-theoretic quantity that averages the negative logs of the
#' probability density functions for \eqn{ŷ} (predicted by the GLM), evaluated
#' at \eqn{y}.
#'
#' Note that cross entropy simply parameterizes prediction error in terms of
#' relative likelihood; if applied to the training set instead of the test set,
#' MCE equals the negative log-likelihood of the model divided by the number of
#' observations. (For logistic regression, this is also known as "log-loss".)
#' Given that GLMs minimize the negative log-likelihood, cross entropy defines a
#' unified loss function for both fitting and cross-validating GLMs.
#'
#' @section R-squared and deviance explained:
#' The prediction R-squared provides a metric analagous to R-squared, except
#' instead of describing how well a model fits the original training data, it
#' describes how well the model predicts independent test data. Instead of using
#' the traditional formula, \eqn{1 - RSS / TSS}, where RSS corresponds to
#' the sum of the squared residuals, \eqn{(y - ŷ)^2}, and TSS is the total sum
#' of squares, it is calculated as \eqn{1 - dev_PRED / dev_NULL}, where
#' \eqn{dev_PRED} is the deviance for the model prediction and \eqn{dev_NULL} is
#' the deviance for the null model. This yields a quantity equivalent to the
#' traditional R-squared formula for linear models with normal error
#' distributions, but a different quantity for exponential family regression
#' models (e.g., logistic regression) that preserves the interpretation of
#' R-squared as the fraction of uncertainty explained (Cameron & Windmeijer,
#' 1997). (See \code{\link{r2d}} for more details.)
#'
#' Note that the prediction R-squared can be negative if overfitting is severe;
#' i.e., when predicting new data, the model performs worse than if one were to
#' always predict the mean response.
#'
#' @section Supported Probability Distributions:
#' Currently the following "families" are supported: Gaussian, binomial,
#' Poisson, negative binomial, zero-inflated Poisson, and zero-inflated negative
#' binomial.
#'
#' @param object A model object of class "glm", "negbin", or "zeroinfl".
#'
#' @param test_data A data frame containing new data for the response and all
#' predictors that were used to fit the model \code{object}.
#'
#' @return A list giving the mean squared error, mean cross entropy, and
#' deviance R-squared between model predictions and actual observations.
#'
#' @seealso \code{\link{r2d}}, \code{\link[stats]{logLik}},
#' \code{\link{deviance.zeroinfl}}
#'
#' @import stats
#' @export
predict_metrics <- function(object, test_data = NULL){
  if(is.null(object$model))
    stop("Model frame missing. Refit model with arg 'model = TRUE'")
  if(is.null(test_data)) test_data <- object$model
  model_type <- class(object)[1]
  family <- NULL
  if(model_type == "lm") family <- "gaussian"
  if(model_type == "glm") family <- object$family$family
  if(model_type == "negbin") family <- "negbin"
  if(model_type == "zeroinfl"){
    family <- if(object$dist == "poisson") "zip" else "zinb"
  }
  if(is.null(family)) stop(paste(model_type, "class not supported"))
  if(model_type == "glm" && !family %in% c("gaussian", "binomial", "poisson"))
    stop(paste(family, "family not supported"))
  y <- unlist(test_data[, names(object$model)[1]])
  if(is.factor(y)) y <- as.integer(y) - 1
  y_hat <- predict(object, test_data, type="response")
  y <- y[!is.na(y_hat)]
  y_hat <- y_hat[!is.na(y_hat)]
  phi <- theta <- NULL
  if(family %in% c("zip", "zinb")){
    y_hat <- predict(object, test_data, type = "count")
    phi <- predict(object, test_data, type = "zero")
  }
  if(family %in% c("negbin","zinb")) theta <- object$theta
  predict_metrics_(y, y_hat, family, phi, theta)
}

predict_metrics_ <- function(y, y_hat, family, phi = NULL, theta = NULL){
  y_bar <- mean(y)
  sigma <- sqrt(mean((y_hat-y)^2))
  N <- length(y)
  ll_null <- sum(
    switch(
      family,
      gaussian = dnorm(y, mean = y_bar, sd = sigma, log = TRUE),
      binomial = dbinom(y, size = 1, prob = y_bar, log = TRUE),
      poisson = dpois(y, lambda = y_bar, log = TRUE),
      negbin = dnbinom(y, size = theta, mu = y_bar, log = TRUE),
      zip = dpois(y, lambda = y_bar, log = TRUE),
      zinb = dnbinom(y, size = theta, mu = y_bar, log = TRUE)
    )
  )
  ll_predicted <- sum(
    switch(
      family,
      gaussian = dnorm(y, mean = y_hat, sd = sigma, log = TRUE),
      binomial = dbinom(y, size = 1, prob = y_hat, log = TRUE),
      poisson = dpois(y, lambda = y_hat, log = TRUE),
      negbin = dnbinom(y, mu = y_hat, size = theta, log = TRUE),
      zip = VGAM::dzipois(y, lambda = y_hat, pstr0 = phi, log = TRUE),
      zinb = VGAM::dzinegbin(y, size = theta, munb = y_hat, pstr0 = phi,
                             log = TRUE)
    )
  )
  ll_saturated <- sum(
    switch(
      family,
      gaussian = -N/2 * log(2*pi*sigma^2),
      binomial = 0,
      poisson = dpois(y, lambda = y, log = TRUE),
      negbin = dnbinom(y, mu = y, size = theta, log = TRUE),
      zip = dpois(y, lambda = y, log = TRUE),
      zinb = dnbinom(y, size = theta, mu = y, log = TRUE)
    )
  )
  dev_pred <- 2 * (ll_saturated - ll_predicted)
  dev_null <- 2 * (ll_saturated - ll_null)
  structure(list(mean_cross_entropy = -ll_predicted / N,
                 mean_squared_error = mean((y - y_hat)^2),
                 R_squared = 1 - dev_pred / dev_null),
            class = "prediction_metrics")
}

