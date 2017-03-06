#' R-squared as Deviance Explained
#'
#' \code{r2d} calculates R-squared as the fraction of deviance explained, which
#' equals traditional R-squared (fraction of variance explained) for linear
#' regression models, but, unlike traditional R-squared, generalizes to
#' exponential family regression models. With optional additional arguments,
#' \code{r2d} can also returns a predictive R-squared for how well the model
#' predicts responses for new observations and/or a cross-validated R-squared
#' with a bootstrapped confidence interval.
#'
#' For standard linear regression models fit with \code{\link[stats]{lm}}, the
#' familiar coefficient of determination, R-squared, can be obtained with
#' \code{\link[stats]{summary.lm}}. However, R-squared is not provided for
#' generalized linear models (GLMs) fit with \code{\link[stats]{glm}},
#' presumably because it can have undesirable properties when applied to GLMs
#' with non-normal error distributions (e.g., binomial, Poisson, etc.). For
#' example, R-squared is no longer guaranteed to lie within the [0,1] interval
#' or to uniformly increase as more predictors are added. Most importantly, the
#' interpretation of R-squared as the fraction of uncertainty explained by the
#' model does not generally hold for exponential family regression models.
#'
#' Cameron and Windmeijer (1997) proposed an R-squared measure (termed
#' \eqn{R_{KL}^2}) for GLMs based on Kullback-Leibler (KL) divergence
#' (\eqn{entropy -  cross-entropy}), which restores all of the desirable
#' properties of R-squared, including its interpretation as the fraction of
#' uncertainty explained. Owing to an equivalence between KL-divergence and
#' deviance, others have referred to this metric as \eqn{R_D^2} (Martin & Hall,
#' 2016), or deviance explained. It is defined as \eqn{1-dev/nulldev}, where
#' \eqn{dev} is the deviance: \eqn{2*(loglik_sat - loglik_fit)}, where
#' \eqn{loglik_sat} is the log-likelihood for the saturated model (a model
#' with 1 free parameter per observation) and \eqn{loglik_fit} is the log-
#' likelihood for the fitted model; and where \eqn{nulldev} is the null
#' deviance: \eqn{2*(loglik_sat - loglik_null)}, where \eqn{loglik_null} is the
#' log-likelihood for the intercept-only model. Following the reasoning of
#' Martin & Hall (2016), the saturated and null models for zero-inflated (ZI)
#' regression are equivalent to the saturated and null models for the
#' corresponding non-ZI regression: e.g., the saturated model for Poisson
#' regression defines the saturated model for ZI-Poisson regression.
#'
#' @param object Fitted model
#'
#' @param newdata An optional data frame in which to look for variables with
#' which to calculate a predictive R-squared
#'
#' @param cv A switch indicating if a predictive R-squared should be estimated
#' by cross-validation via a call to \code{\link{cv_r2}}.
#'
#' @param ... Additional arguments to be passed to \code{\link{cv_r2}}.
#'
#' @return Object of class "R2" with the following items:
#' \enumerate{
#'  \item\describe{
#'   \item{R2fit}{an R-squared statistic for the model fit}}
#'  \item\describe{
#'   \item{R2new}{if \code{newdata} was provided, an R-squared
#'   statistic for how well the model predicts new data}}
#'  \item\describe{
#'   \item{R2cv}{if \code{cv = TRUE}, an object returned by \code{\link{cv_r2}}}
#'   }}
#'
#' @seealso \code{\link{predict_metrics}}, \code{\link{deviance.zeroinfl}},
#' \code{\link{cv_r2}}
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
r2d <- function (object, newdata = NULL, cv = FALSE, ...) {
  model_type <- class(object)[1]
  R2fit <- switch(model_type,
                  lm = 1 - var(resid(object)) / var(object$model[[1]]),
                  glm = 1 - object$deviance / object$null.deviance,
                  negbin = 1 - object$deviance / object$null.deviance,
                  zeroinfl = predict_metrics(object)$R_squared)
  R2new <- NULL
  if(!is.null(newdata)){
    R2new <- predict_metrics(object, newdata)$R_squared
  }
  R2cv <- NULL
  if(cv){
    R2cv <- cv_r2(object, ...)
  }
  structure(list(R2fit = R2fit, R2new = R2new, R2cv = R2cv), class = "R2")
}
