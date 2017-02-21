#' Proportion of variance and deviance explained
#'
#' Returns the fraction of variance explained and the fraction of deviance
#' explained by a fitted model object.
#'
#' For standard linear regression models fit with \code{\link[stats]{lm}}, the
#' familiar coefficient of determination, R-squared, can be obtained with
#' \code{\link[stats]{summary.lm}}. However, R-squared is not provided for
#' generalized linear models (GLMs) fit with \code{\link[stats]{glm}},
#' presumably because it can have undesirable properties when applied to GLMs
#' if the error distribution is non-normal (e.g., binomial, Poisson, etc.).
#' For example, it can lie outside the [0,1] interval and decrease as predictors
#' are added. Most importantly, the interpretation of R-squared as the fraction
#' of uncertainty explained by the model does not generally hold for exponential
#' family regression models.
#'
#' Cameron and Windmeijer (1997) proposed an R-squared measure for GLMs based on
#' Kullback-Leibler (KL) divergence (\eqn{entropy -  cross-entropy}), which they
#' termed \eqn{R_{KL}^2}, that restores all of the desirable properties of
#' R-squared, including its interpretation as the fraction of uncertainty
#' explained. Owing to an equivalence between KL-divergence and deviance, others
#' have referred to as \eqn{R_D^2} (Martin & Hall, 2016), \eqn{D^2}
#' (\code{\link[modEvA]{Dsquared}}), or \code{dev.ratio}
#' (\code{\link[glmnet]{glmnet}}). (\code{r2d2} adopts the latter terminology,
#' and the \code{\link{summary.beset_glm}} and
#' \code{\link{summary.beset_zeroinfl}} methods refer to this measure as
#' "Deviance Explained".
#'
#' The \code{dev_ratio} is defined to be \eqn{1-dev/nulldev}, where \eqn{dev} is
#' the deviance: \eqn{2*(loglik_sat - loglik_fit)}, where \eqn{loglik_sat} is
#' the log-likelihood for the saturated model (a model witha free parameter per
#' observation) and \eqn{loglik_fit} is the log-likelihood for the fitted model;
#' and where \eqn{nulldev} is the null deviance:
#' \eqn{2*(loglik_sat - loglik_null)}, where \eqn{loglik_null} is the
#' log-likelihood for the intercept-only model.
#'
#' @return \code{r2d2} returns a list with both variance explained (named
#' \code{"r_squared"}) and deviance explained (named \code{"dev_ratio"}) for
#' GLMs and zero-inflated models. For GLMs with a normal error distribution, the
#' two are equivalent. For all other GLMs, the deviance explained provides
#' the more meaningful measure of model fit.
#'
#' @references
#' Colin Cameron A, Windmeijer FAG. (1997) An R-squared measure of goodness
#' of fit for some common nonlinear regression models. J Econometrics, 77,
#' 329–342.
#'
#' Jacob Martin & Daniel B. Hall (2016) R2 measures for zero-inflated regression
#' models for count data with excess zeros, Journal of Statistical Computation
#' and Simulation, 86:18, 3777-3790, DOI: 10.1080/00949655.2016.1186166
#'
#' @export
r2d2 <- function (x, ...) {
  UseMethod("r2d2", x)
}
#' @export
r2d2.lm <- function(object){
  r2 <- 1 - var(resid(object, type = "response")) / var(object$model[[1]])
  structure(list(r_squared = r2, dev_ratio = r2), class = "r2d2")
}
#' @export
r2d2.glm <- function(object){
  y <- object$model[[1]]
  if(is.factor(y)) y <- as.integer(y) - 1
  r2 <- 1 - var(resid(object, type = "response")) / var(y)
  d2 <- 1 - object$deviance / object$null.deviance
  structure(list(r_squared = r2, dev_ratio = d2), class = "r2d2")
}
#' @export
r2d2.zeroinfl <- function(object){
  r2 <- 1 - var(resid(object, type = "response")) / var(object$model[[1]])
  null <- glm(object$model[[1]] ~ 1, family = "poisson")
  d2 <- 1 - deviance(object) / deviance(null)
  structure(list(r_squared = r2, dev_ratio = d2), class = "r2d2")
}

#' @export
r2d2.beset_glm <- function(object){
  r2d2(object$best_model)
}
