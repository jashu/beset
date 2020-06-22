#' Prediction Metrics
#'
#' Calculate deviance, mean absolute error, mean cross entropy, mean squared
#' error, and R-squared (as fraction of deviance explained) to measure how well
#' a generalized linear model predicts test data.
#'
#' \code{predict_metrics} uses a generalized linear model previously fit to a
#' training set to predict responses for a test set. It then computes the
#' residual sum of squares and the log-likelihood of the predicted responses, as
#' well as the log-likelihoods for the corresponding saturated model (a model
#' with one free parameter per observation) and null model (a model with only an
#' intercept) fit to the true responses. These quantities are used to derive the
#' metrics described below.
#'
#' \code{predict_metrics_} is the workhorse function: there is no need to call
#' this directly but is more efficient for programming if the actual response
#' vector, the predicted response vector(s), and dispersion parameter (if
#' applicable) have already been calculated.
#'
#' @section Mean squared error (MSE) and mean cross entropy (MCE):
#' MSE is the average of the squared difference between each predicted response
#' \eqn{ŷ} and the actual observed response \eqn{y}. MCE is an analagous
#' information-theoretic quantity that averages the negative logs of the
#' probability density functions for \eqn{ŷ} evaluated at \eqn{y}.
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
#' @param y Vector of actual observed values
#'
#' @param y_hat Vector of predicted values.
#'
#' @param family Character string giving the name of the error distribution for
#'   \code{y_hat}. Currently supported distributions are "gaussian", "binomial",
#'   "poisson", "negbin" (negative binomial), "zip" (zero-inflated Poisson), and
#'    "zinb" (zero-inflated negative binomial).
#'
#' @param mu (Required for zero-inflated models) the predicted values for the
#' count component of the model.
#'
#' @param phi (Required for zero-inflated models) the predicted values for the
#'  zero component of the model.
#'
#' @param theta (Required for negative binomial models) the estimated dispersion
#' parameter \code{theta}.
#'
#' @return For binomial models, a list giving the area under the ROC curve
#' ("auc"), mean cross entropy AKA log loss ("mce"), mean squared error AKA
#' Brier score ("mse"), and deviance R-squared ("rsq'). For all other models,
#' mean absolute error ("mae") takes the place of "auc" but otherwise the
#' metrics are the same.
#'
#' @seealso \code{\link{r2d}}, \code{\link[stats]{logLik}}
#'
#' @importFrom stats binomial dnbinom dpois pnbinom pnorm poisson ppois predict
#'  qnorm
#' @export
predict_metrics <- function(object, test_data){
  model_type <- class(object)[1]
  family <- switch(model_type,
                   lm = "gaussian",
                   glm = object$family$family,
                   negbin = "negbin",
                   zeroinfl = if(object$dist == "poisson") "zip" else "zinb")
  if(is.null(family)) stop(paste(model_type, "class not supported"))
  if(model_type == "glm" && !family %in% c("gaussian", "binomial", "poisson"))
    stop(paste(family, "family not supported"))
  y <- if(model_type == "zeroinfl"){
    purrr::as_vector(test_data[[all.vars(object$terms$full)[1]]])
  } else {
    purrr::as_vector(test_data[[all.vars(object$terms)[1]]])
  }
  if(is.null(y)) stop("Response variable not found in `test_data`")
  y_hat <- stats::predict(object, test_data, type="response")
  mu <- phi <- theta <- NULL
  if(family %in% c("zip", "zinb")){
    mu <- stats::predict(object, test_data, type = "count")
    phi <- stats::predict(object, test_data, type = "zero")
  }
  if(family %in% c("negbin","zinb")) theta <- object$theta
  predict_metrics_(y, y_hat, family, mu, phi, theta)
}

#' @rdname predict_metrics
#' @export
predict_metrics_ <- function(
  y, y_hat, family, theta = NULL, mu = NULL, phi = NULL
){
  na <- which(is.na(y) | is.na(y_hat))
  if(length(na) > 0){
    y <- y[-na]
    y_hat <- y_hat[-na]
    mu <- mu[-na]
    phi <- phi[-na]
  }
  if(is.factor(y)) y <- as.integer(y) - 1
  y_bar <- mean(y)
  sigma <- sqrt(mean((y_hat-y)^2))
  N <- length(y)
  ll_null <- sum(
    switch(
      family,
      gaussian = stats::dnorm(y, mean = y_bar, sd = sigma, log = TRUE),
      binomial = stats::dbinom(y, size = 1, prob = y_bar, log = TRUE),
      poisson = stats::dpois(y, lambda = y_bar, log = TRUE),
      negbin = stats::dnbinom(y, size = theta, mu = y_bar, log = TRUE),
      zip = stats::dpois(y, lambda = y_bar, log = TRUE),
      zinb = stats::dnbinom(y, size = theta, mu = y_bar, log = TRUE)
    )
  )
  ll_predicted <- sum(
    switch(
      family,
      gaussian = stats::dnorm(y, mean = y_hat, sd = sigma, log = TRUE),
      binomial = stats::dbinom(y, size = 1, prob = y_hat, log = TRUE),
      poisson = stats::dpois(y, lambda = y_hat, log = TRUE),
      negbin = stats::dnbinom(y, mu = y_hat, size = theta, log = TRUE),
      zip = VGAM::dzipois(y, lambda = mu, pstr0 = phi, log = TRUE),
      zinb = VGAM::dzinegbin(y, size = theta, munb = mu, pstr0 = phi,
                             log = TRUE)
    )
  )
  ll_saturated <- sum(
    switch(
      family,
      gaussian = -N/2 * log(2*pi*sigma^2),
      binomial = 0,
      poisson = stats::dpois(y, lambda = y, log = TRUE),
      negbin = stats::dnbinom(y, size = theta, mu = y,  log = TRUE),
      zip = stats::dpois(y, lambda = y, log = TRUE),
      zinb = stats::dnbinom(y, size = theta, mu = y, log = TRUE)
    )
  )
  dev_pred <- 2 * (ll_saturated - ll_predicted)
  if(family == "gaussian") dev_pred <- dev_pred * sigma^2
  dev_null <- 2 * (ll_saturated - ll_null)
  if(family == "gaussian") dev_null <- dev_null * sigma^2
  auc <- NULL
  if(family == "binomial"){
    x1 <- y_hat[y==1]; n1 <- length(x1);
    x2 <- y_hat[y==0]; n2 <- length(x2);
    r <- rank(c(x1, x2))
    auc <- if(n1 > 0 && n2 > 0){
      (sum(r[1:n1]) - n1 * (n1 + 1) / 2) / n1 / n2
        } else NA_real_
  }
  structure(
    c(if(!is.null(auc)) list(auc = auc) else list(mae = mean(abs(y_hat - y))),
    list(mce = -ll_predicted / N, mse = sigma^2,
         rsq = if(length(y) > 2) 1 - dev_pred/dev_null else NA_real_)),
    class = "prediction_metrics", family = family, theta = theta
  )
}
