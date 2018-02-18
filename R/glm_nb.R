#' Fit a Negative Binomial Generalized Linear Model
#'
#' A modification of \code{\link[MASS]{glm.nb}} that can be programmed with.
#' The main change is to avoid non-standard evaluation of the link argument so
#' that this argument is passed as a "quoted" string, which requires the
#' removal of the line `link <- substitute(link)`. This in turn enables the
#' removal of `do.call` syntax. Other edits were done to define the namespaces
#' of other MASS functions that are called by `glm.nb`, to change the default
#' arguments to avoid returning copies of the model data and the response
#' vector, and to return \code{$data} to maintain consistency with the object
#' returned by \code{\link[stats]{glm}}.
#'
#' @inheritParams stats::glm
#'
#' @param init.theta Optional initial value for the theta parameter. If omitted
#' a moment estimator after an initial fit using a Poisson GLM is used.
#'
#' @param link The link function. Currently must be one of "log", "sqrt" or
#' "identity".
#'
#' @export

glm_nb <- function (formula, data, weights, subset, na.action, start = NULL,
                    etastart, mustart, control = stats::glm.control(...),
                    method = "glm.fit", model = FALSE, x = FALSE, y = FALSE,
                    contrasts = NULL, ..., init.theta, link = "log"){
  loglik <- function(n, th, mu, y, w){
    sum(w * (lgamma(th + y) - lgamma(th) - lgamma(y + 1) + th * log(th) + y *
               log(mu + (y == 0)) - (th + y) * log(th + mu)))
  }

  fam0 <- if (missing(init.theta)){
    stats::poisson(link = link)
  } else  {
    MASS::negative.binomial(theta = init.theta, link = link)
  }
  mf <- Call <- match.call()
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "etastart", "mustart", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval.parent(mf)
  Terms <- attr(mf, "terms")
  if (method == "model.frame") return(mf)
  Y <- stats::model.response(mf, "numeric")
  X <- if (!stats::is.empty.model(Terms)){
    stats::model.matrix(Terms, mf, contrasts)
  } else {
    matrix(nrow = NROW(Y), ncol = 0)
  }
  w <- stats::model.weights(mf)
  if (!length(w)){
    w <- rep(1, nrow(mf))
  } else if (any(w < 0)) stop("negative weights not allowed")
  offset <- stats::model.offset(mf)
  mustart <- stats::model.extract(mf, "mustart")
  etastart <- stats::model.extract(mf, "etastart")
  n <- length(Y)
  if (!missing(method)) {
    if (!exists(method, mode = "function"))
      stop(gettextf("unimplemented method: %s", sQuote(method)),
           domain = NA)
    glm.fitter <- get(method)
  } else {
    method <- "glm.fit"
    glm.fitter <- stats::glm.fit
  }
  if (control$trace > 1) message("Initial fit:")
  fit <- glm.fitter(x = X, y = Y, w = w, start = start, etastart = etastart,
                    mustart = mustart, offset = offset, family = fam0,
                    control = list(maxit = control$maxit,
                                   epsilon = control$epsilon,
                                   trace = control$trace > 1),
                    intercept = attr(Terms, "intercept") > 0)
  class(fit) <- c("glm", "lm")
  mu <- fit$fitted.values
  th <- as.vector(MASS::theta.ml(Y, mu, sum(w), w, limit = control$maxit,
                                 trace = control$trace > 2))
  if (control$trace > 1){
    message(gettextf("Initial value for 'theta': %f", signif(th)), domain = NA)
  }
  fam <- MASS::negative.binomial(theta = th, link = link)
  iter <- 0
  d1 <- sqrt(2 * max(1, fit$df.residual))
  d2 <- del <- 1
  g <- fam$linkfun
  Lm <- loglik(n, th, mu, Y, w)
  Lm0 <- Lm + 2 * d1
  while ((iter <- iter + 1) <= control$maxit &&
         (abs(Lm0 - Lm)/d1 + abs(del)/d2) > control$epsilon) {
    eta <- g(mu)
    fit <- glm.fitter(x = X, y = Y, w = w, etastart = eta, offset = offset,
                      family = fam, control = list(maxit = control$maxit,
                                                   epsilon = control$epsilon,
                                                   trace = control$trace > 1),
                      intercept = attr(Terms, "intercept") > 0)
    t0 <- th
    th <- MASS::theta.ml(Y, mu, sum(w), w, limit = control$maxit,
                         trace = control$trace > 2)
    fam <- MASS::negative.binomial(theta = th, link = link)
    mu <- fit$fitted.values
    del <- t0 - th
    Lm0 <- Lm
    Lm <- loglik(n, th, mu, Y, w)
    if (control$trace) {
      Ls <- loglik(n, th, Y, Y, w)
      Dev <- 2 * (Ls - Lm)
      message(sprintf("Theta(%d) = %f, 2(Ls - Lm) = %f",
                      iter, signif(th), signif(Dev)), domain = NA)
    }
  }
  if (!is.null(attr(th, "warn")))
    fit$th.warn <- attr(th, "warn")
  if (iter > control$maxit) {
    warning("alternation limit reached")
    fit$th.warn <- gettext("alternation limit reached")
  }
  if (length(offset) && attr(Terms, "intercept")) {
    null.deviance <- if (length(Terms))
      glm.fitter(X[, "(Intercept)", drop = FALSE], Y, w, offset = offset,
                 family = fam, control = list(maxit = control$maxit,
                                              epsilon = control$epsilon,
                                              trace = control$trace > 1),
                 intercept = TRUE)$deviance
    else fit$deviance
    fit$null.deviance <- null.deviance
  }
  class(fit) <- c("negbin", "glm", "lm")
  fit$terms <- Terms
  fit$formula <- as.vector(attr(Terms, "formula"))
  Call$init.theta <- signif(as.vector(th), 10)
  Call$link <- link
  fit$call <- Call
  if (model)
    fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  if (x) fit$x <- X
  if (!y) fit$y <- NULL
  fit$theta <- as.vector(th)
  fit$SE.theta <- attr(th, "SE")
  fit$twologlik <- as.vector(2 * Lm)
  fit$aic <- -fit$twologlik + 2 * fit$rank + 2
  fit$contrasts <- attr(X, "contrasts")
  fit$xlevels <- stats::.getXlevels(Terms, mf)
  fit$method <- method
  fit$control <- control
  fit$offset <- offset
  fit$data <- data
  fit
}

logLik.negbin <- function (object, ...) {
  if (nargs() > 1L)
    warning("extra arguments discarded")
  p <- object$rank + 1L
  val <- object$twologlik/2
  attr(val, "df") <- p
  attr(val, "nobs") <- sum(!is.na(object$residuals))
  class(val) <- "logLik"
  val
}
