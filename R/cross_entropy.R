#' Cross Entropy
#'
#' Calculate mean cross entropy to measure how well a generalized linear model
#' predicts test data.
#'
#' \code{cross_entropy} uses a generalized linear model (GLM) previously fit to
#' a training set to predict responses for a test set. It then computes the
#' cross entropy between each predicted response \eqn{ŷ} and the corresponding
#' observed response \eqn{y} as the negative log of the probability density
#' function for \eqn{ŷ} (predicted by the GLM) evaluated at \eqn{y}, and returns
#' the average cross entropy over all test examples.
#'
#' Note that cross entropy simply parameterizes prediction error in terms of
#' relative likelihood; if applied to the training set instead of the test set,
#' the value returned by \code{cross_entropy} equals the negative log-likelihood
#' of the model divided by the number of observations. (For logistic regression,
#' this is also known as "log-loss".) Given that GLMs minimize the negative
#' log-likelihood, \code{cross_entropy} defines a unified loss function for both
#' fitting and cross-validating GLMs.
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
cross_entropy <- function(object, test_data = NULL){
  model_type <- class(object)[1]
  if(!model_type %in% c("glm", "negbin", "zeroinfl"))
    stop(paste(model_type, "class not supported"))
  if(model_type == "glm"){
    if(!object$family$family %in% c("gaussian", "binomial", "poisson"))
      stop(paste(object$family$family, "family not supported"))
  }
  if(is.null(test_data)) test_data <- object$model
  y_obs <- unlist(test_data[, names(object$model)[1]])
  y_hat <- predict(object, test_data, type="response")
  if(class(object) == "zeroinfl"){
    if(object$dist == "poisson") family <- "zip" else family <- "zinb"
  } else family <- object$family$family
  if(grepl("Negative Binomial", family)) family <- "negbin"
  log_prob <- switch(family,
                 gaussian = dnorm(y_obs, mean = y_hat,
                                  sd = sigma(object), log = TRUE),
                 binomial = y_obs * log(y_hat) + (1 - y_obs) * log(1 - y_hat),
                 poisson = dpois(y_obs, lambda = y_hat, log = TRUE),
                 negbin = dnbinom(y_obs, mu = y_hat, size=object$theta,
                                  log = TRUE),
                 zip = VGAM::dzipois(y_obs,
                                     lambda = predict(object,
                                                      test_data,
                                                      type = "count"),
                                     pstr0 = predict(object,
                                                     test_data,
                                                     type = "zero"),
                                     log = TRUE),
                 zinb = VGAM::dzinegbin(y_obs,
                                        size = object$theta,
                                        munb = predict(object,
                                                       test_data,
                                                       type = "count"),
                                        pstr0 = predict(object,
                                                        test_data,
                                                        type = "zero"),
                                        log = TRUE))
  -1 * mean(log_prob)
}
