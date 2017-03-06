#' Deviance Method for Zero-Inflated Models
#'
#' Implements a \code{deviance} method for the \code{"zeroinfl"} class.
#'
#' The \code{pscl} package does not provide a method for calculating the
#' deviance of zero-inflated (ZI) models. This method is added here with
#' deviance calculated in the ususal way: eqn{2*(loglik_sat - loglik_fit)},
#' where \eqn{loglik_sat} is the log-likelihood for the saturated model and
#' \eqn{loglik_fit} is the log-likelihood for the fitted model. The \code{pscl}
#' package defines the log-likelihood for fitted ZI models. Following the
#' reasoning of Martin & Hall (2016), the saturated model for ZI regression is
#' taken to be equivalent to the saturated model for the corresponding non-ZI
#' regression: e.g., the saturated model for ZI-Poisson regression is the same
#' as the saturated model for Poisson regression.
#'
#' @inheritParams stats::deviance
#'
#' @return The value of the deviance extracted from the model \code{object}.
#'
#' @references Jacob Martin & Daniel B. Hall (2016) R2 measures for
#' zero-inflated regression models for count data with excess zeros, Journal of
#' Statistical Computation and Simulation, 86:18, 3777-3790,
#' DOI: 10.1080/00949655.2016.1186166
#'
#' @export
deviance.zeroinfl <- function(object, ...){
  loglik_sat <- sum(
    switch(
      object$dist,
      poisson = dpois(object$y, lambda = object$y, log = TRUE),
      negbin = dnbinom(object$y, mu = object$y, size=object$theta, log = TRUE)
      )
    )
  as.numeric(2 * (loglik_sat - logLik(object)))
}
