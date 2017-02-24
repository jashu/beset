#' Prediction Metrics
#'
#' Calculate squared-residual (mean squared error and R-squared) and analagous
#' deviance-residual (mean cross entropy and deviance explained) metrics to
#' measure how well a generalized linear model predicts test data.
#'
#' \code{prediction_metrics} uses a generalized linear model (GLM) previously
#' fit to a training set to predict responses for a test set. It then computes
#' the residual sum of squares and the log-likelihood of the predicted
#' responses, as well as the log-likelihoods for the corresponding saturated
#' model (a model with one free parameter per observation) and null model (a
#' model with only an intercept) fit to the true responses. These quantities
#' are used to derive the metrics described in the following sections.
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
#' describes how well the model predicts independent test data. It is calculated
#' using the traditional formula \eqn{1 - RSS / TSS}, where RSS corresponds to
#' the sum of the squared residuals (\eqn{(y - ŷ)^2}) and TSS is the total sum
#' of squares.
#'
#' The traditional interpretation of R-squared as the fraction of uncertainty
#' explained by the model does not generally hold for exponential family
#' regression models. Deviance explained, originally termed \eqn{R_{KL}^2} by
#' Cameron and Windmeijer (1997), is an R-squared-analagous measure based
#' on Kullback-Leibler (KL) divergence (\eqn{entropy -  cross-entropy}) that
#' does retain this interpretation for the wider class of GLMs. Owing to an
#' equivalence between KL-divergence and deviance, others have referred to
#' this metric as \eqn{R_D^2} (Martin & Hall, 2016), \eqn{D^2}
#' (\code{\link[modEvA]{Dsquared}}), or \code{dev.ratio}
#' (\code{\link[glmnet]{glmnet}}). It is calculated as
#' \eqn{1 - dev_PRED / dev_NULL}, where \eqn{dev_PRED} is the deviance for the
#' model prediction and \eqn{dev_NULL} is the deviance for the null model.

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
#' @return A number giving the cross entropy error between model predictions and
#' actual observations.
#'
#' @seealso \code{\link{beset_glm}}, \code{\link{beset_zeroinfl}}
#'
#' @export
prediction_metrics <- function(object, test_data = NULL){
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
  y_bar <- mean(y)
  N <- length(y)
  if(family == "gaussian"){
    n <- nrow(object$model)
    sd <- sigma(object) * sqrt((n-dim(model.matrix(object))[2])/n)
  }
  ll_null <- sum(
    switch(
      family,
      gaussian = dnorm(y, mean = y_bar, sd = sd, log = TRUE),
      binomial = dbinom(y, size = 1, prob = y_bar, log = TRUE),
      poisson = dpois(y, lambda = y_bar, log = TRUE),
      negbin = dnbinom(y, size=object$theta, mu = y_bar, log = TRUE),
      zip = dpois(y, lambda = y_bar, log = TRUE),
      zinb = dnbinom(y, size=object$theta, mu = y_bar, log = TRUE)
    )
  )
  ll_predicted <- sum(
    switch(
      family,
      gaussian = dnorm(y, mean = y_hat, sd = sd, log = TRUE),
      binomial = dbinom(y, size = 1, prob = y_hat, log = TRUE),
      poisson = dpois(y, lambda = y_hat, log = TRUE),
      negbin = dnbinom(y, mu = y_hat, size=object$theta, log = TRUE),
      zip = VGAM::dzipois(y,
                          lambda = predict(object, test_data, type = "count"),
                          pstr0 = predict(object, test_data, type = "zero"),
                          log = TRUE),
      zinb = VGAM::dzinegbin(y, size = object$theta,
                             munb = predict(object, test_data, type = "count"),
                             pstr0 = predict(object, test_data, type = "zero"),
                             log = TRUE)
    )
  )
  ll_saturated <- sum(
    switch(
      family,
      gaussian = -length(y)/2 * log(2*pi*sd^2),
      binomial = 0,
      poisson = dpois(y, lambda = y, log = TRUE),
      negbin = dnbinom(y, mu = y, size=object$theta, log = TRUE),
      zip = dpois(y, lambda = y, log = TRUE),
      zinb = dnbinom(y, size=object$theta, mu = y, log = TRUE)
    )
  )
  KLd_pred <- 2 * (ll_saturated - ll_predicted)
  KLd_null <- 2 * (ll_saturated - ll_null)
  structure(list(mean_cross_entropy = -ll_predicted / N,
                 mean_squared_error = mean((y - y_hat)^2),
                 R_squared = 1 - KLd_pred / KLd_null),
            class = "prediction_metrics")
}
