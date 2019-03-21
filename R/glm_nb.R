#' Fit a Negative Binomial Generalized Linear Model
#'
#' A modification of \code{\link[MASS]{glm.nb}} to provide a more efficient
#' workhorse function analagous to \code{\link[stats]{glm.fit}} where the
#' response vector, design matrix, and family have already been calculated.
#'
#' @inheritParams beset_glm
#' @param x A design matrix of dimension \code{n * p}
#' @param y A vector of observations of length \code{n}.
#' @param offset A vector of length \code{nobs} specifying an \emph{a priori}
#'               known component that will be added to the linear predictor
#'               before applying the link function. Default is NULL.
#' @param family The result of a call to either \code{\link[stats]{poisson}} or
#'               \code{\link[MASS]{negative.binomial}}.
#' @param control A list of parameters for controlling the fitting process to be
#'                passed to \code{\link[stats]{glm.control}}
#' @export

glm_nb <- function(x, y, weights = rep(1, nobs), start = NULL, etastart = NULL,
                   mustart = NULL, offset = rep(0, nobs), family = poisson(),
                   control = list(), intercept = TRUE){
  loglik <- function(th, mu, y, w){
    sum(w * (lgamma(th + y) - lgamma(th) - lgamma(y + 1) + th * log(th) + y *
               log(mu + (y == 0)) - (th + y) * log(th + mu)))
  }
  if(!grepl("(poisson)|(Negative Binomial)", family$family)){
    stop(paste(
      "`family` argument must be an object returned by",
      "`poisson()` or `negative.binomial()`", sep = "\n")
    )
  }
  control <- do.call("glm.control", control)
  xnames <- dimnames(x)[[2L]]
  ynames <- if (is.matrix(y)) rownames(y) else names(y)
  conv <- FALSE
  nobs <- NROW(y)
  nvars <- ncol(x)
  fam0 <- family
  w <- if(is.null(weights)) rep.int(1, nobs) else weights
  if(is.null(offset)) offset <- rep.int(0, nobs)
  fit <- stats::glm.fit(x = x, y = y, w = w, start = start, etastart = etastart,
                    mustart = mustart, offset = offset, family = fam0,
                    control = control, intercept = intercept)
  class(fit) <- c("glm", "lm")
  mu <- fit$fitted.values
  th <- as.vector(MASS::theta.ml(y, mu, sum(w), w, limit = control$maxit))
  fam <- MASS::negative.binomial(theta = th, link = family$link)
  iter <- 0
  d1 <- sqrt(2 * max(1, fit$df.residual))
  d2 <- del <- 1
  g <- fam$linkfun
  Lm <- loglik(th, mu, y, w)
  Lm0 <- Lm + 2 * d1
  while((iter <- iter + 1) <= control$maxit &&
        (abs(Lm0 - Lm)/d1 + abs(del)/d2) > control$epsilon) {
    eta <- g(mu)
    fit <- stats::glm.fit(
      x = x, y = y, w = w, etastart = eta, offset = offset, family = fam,
      control = control, intercept = intercept)
    t0 <- th
    th <- MASS::theta.ml(y, mu, sum(w), w, limit = control$maxit)
    fam <- MASS::negative.binomial(theta = th, link = family$link)
    mu <- fit$fitted.values
    del <- t0 - th
    Lm0 <- Lm
    Lm <- loglik(th, mu, y, w)
  }
  if (!is.null(attr(th, "warn")))
    fit$th.warn <- attr(th, "warn")
  if (iter > control$maxit) {
    warning("alternation limit reached")
    fit$th.warn <- gettext("alternation limit reached")
  }
  if (length(offset) && intercept) {
    fit$null.deviance <- stats::glm.fit(
      x[, "(Intercept)", drop = FALSE], y, w, offset = offset, family = fam,
      control = control, intercept = TRUE)$deviance
  }
  class(fit) <- c("negbin", "glm", "lm")
  fit$theta <- as.vector(th)
  fit$SE.theta <- attr(th, "SE")
  fit$twologlik <- as.vector(2 * Lm)
  fit$aic <- -fit$twologlik + 2 * fit$rank + 2
  fit$contrasts <- attr(x, "contrasts")
  fit$control <- control
  fit$offset <- offset
  fit
}

logLik.negbin <- function (object, ...) {
  p <- object$rank + 1L
  val <- object$twologlik/2
  attr(val, "df") <- p
  attr(val, "nobs") <- sum(!is.na(object$residuals))
  class(val) <- "logLik"
  val
}

theta_ml <- function(y, mu, n = sum(weights), weights, limit = 10,
                     eps = .Machine$double.eps^0.25, trace = FALSE)
{
  score <- function(n, th, mu, y, w){
    sum(
      w * (
        digamma(th + y) - digamma(th) + log(th) + 1 - log(th + mu) -
          (y + th)/(mu + th)
        )
      )
  }
  info <- function(n, th, mu, y, w){
    sum(
      w * (
        -trigamma(th + y) + trigamma(th) - 1/th + 2/(mu + th) -
          (y + th)/(mu + th)^2
        )
      )
  }
  if (inherits(y, "lm")) {
    mu <- y$fitted.values
    y <- if (is.null(y$y))
      mu + residuals(y)
    else y$y
  }
  if (missing(weights))
    weights <- rep(1, length(y))
  t0 <- n/sum(weights * (y/mu - 1)^2)
  it <- 0
  del <- 1
  if (trace)
    message(sprintf("theta.ml: iter %d 'theta = %f'", it,
                    signif(t0)), domain = NA)
  while ((it <- it + 1) < limit && abs(del) > eps) {
    t0 <- abs(t0)
    del <- score(n, t0, mu, y, weights)/(i <- info(n, t0,
                                                   mu, y, weights))
    t0 <- t0 + del
    if (trace)
      message("theta.ml: iter", it, " theta =", signif(t0))
  }
  if (t0 < 0) {
    t0 <- 0
    warning("estimate truncated at zero")
    attr(t0, "warn") <- gettext("estimate truncated at zero")
  }
  if (it == limit) {
    warning("iteration limit reached")
    attr(t0, "warn") <- gettext("iteration limit reached")
  }
  attr(t0, "SE") <- sqrt(1/i)
  t0
}
