#' Predict Methods for \code{beset} Objects
#'
#' @param object A \code{beset_elnet} object.
#'
#' @param type Type of prediction required. Type "link" gives the linear
#' predictors for "binomial" and "poisson" models; for "gaussian" models it
#' gives the fitted values. Type "response" gives the fitted probabilities for
#' "binomial", fitted mean for "poisson", and fitted values for  "gaussian".
#' Type "coefficients" computes the coefficients. Note that for "binomial"
#' models, results are returned only for the class corresponding to the second
#' level of the factor response. Type "class" applies only to "binomial" models,
#' and produces the class label corresponding to the maximum probability. Type
#' "nonzero" returns a list of the indices of the nonzero coefficients.
#'
#' @inheritParams stats::predict.lm
#' @inheritParams glmnet::predict.glmnet
#' @inheritParams get_best

#' @export
predict.beset <- function(object, newdata, type = "response",
                          newoffset = NULL, alpha = NULL, lambda = NULL,
                          n_pred = NULL, metric = "auto", oneSE = TRUE,
                          na.action = na.pass, tt = NULL, ...){
  if(inherits(object, "rf")){
    return(predict.beset_rf(object, newdata, type = "response", ...))
  }
  metric <- tryCatch(
    match.arg(metric, c("auto", "auc", "mae", "mce", "mse", "rsq")),
    error = function(c){
      c$message <- gsub("arg", "metric", c$message)
      c$call <- NULL
      stop(c)
    }
  )
  tryCatch(
    if(
      (metric == "auc" && object$family != "binomial") ||
      (metric == "mae" && object$family == "binomial")
    ) error = function(c){
      c$message <- paste(metric, "not available for", object$family, "models")
      c$call <- NULL
      stop(c)
    }
  )
  if(metric == "auto"){
    metric <- if(object$family == "gaussian") "mse" else "mce"
  }
  if(is.null(tt)) tt <- terms(object)
  if (missing(newdata) || is.null(newdata)) {
    X <- model.matrix(object)
    newoffset <- object$parameters$offset
  }
  else {
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, na.action = na.action,
                     xlev = object$xlevels)
    X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    if("(Intercept)" %in% colnames(X)) X <- X[, -1, drop = FALSE]
    if(is.null(newoffset) && all(object$parameters$fit$offset == 0))
      newoffset <- rep(0, nrow(X))
  }
  if(inherits(object, "elnet")){
    model <- get_best.beset_elnet(object, alpha = alpha, lambda = lambda,
                                  metric = metric, oneSE = oneSE, ...)
    yhat <- glmnet::predict.glmnet(model, newx = X, s = model$best_lambda,
                                   type = type, newoffset = newoffset, ...)
  } else {
    model <- get_best.beset_glm(object, alpha = alpha, lambda = lambda,
                                  metric = metric, oneSE = oneSE, ...)
    yhat <- object$family$linkinv(
      X[, names(coef(object))] %*% coef(object) + newoffset)
  }
  as.vector(yhat)
}

#' @export
model.matrix.beset <- function(object, ...){
  object$parameters$x
}

#' @export
predict.nested <- function(object, newdata, type = "response",
                           newoffset = NULL, alpha = NULL, lambda = NULL,
                           n_pred = NULL, metric = "auto", oneSE = TRUE,
                           na.action = na.pass, ...){
  map(object$beset,
           ~ predict(., newdata, type, newoffset, alpha, lambda, n_pred, metric,
                     oneSE, na.action, object$terms)) %>%
    transpose %>% simplify_all %>% map_dbl(mean)

}

#' @export
predict.beset_rf <- function(object, newdata, type = "response", ...){
  map(object$forests, ~ predict(., newdata, type, ...)) %>%
    transpose %>% simplify_all %>% map_dbl(mean)
}


