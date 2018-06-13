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
predict.beset_elnet <- function(object, newdata, type = "response",
                                newoffset = NULL, alpha = NULL, lambda = NULL,
                                metric = c("mce", "auc", "mae", "mce", "mse"),
                                error = c("auto", "btwn_fold_se",
                                          "btwn_rep_range", "none"),
                                na.action = na.pass, ...){
  tt <- terms(object)
  if (missing(newdata) || is.null(newdata)) {
    X <- model.matrix(object)
    newoffset <- object$glmnet_parameters$offset
  }
  else {
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, na.action = na.action,
                     xlev = object$xlevels)
    X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    if("(Intercept)" %in% colnames(X)) X <- X[,-1]
    if(is.null(newoffset) && all(object$parameters$fit$offset == 0))
      newoffset <- rep(0, nrow(X))
  }
  model <- get_best(object, alpha = alpha, lambda = lambda, metric = metric,
                    error = error, ...)
  yhat <- glmnet::predict.glmnet(model, newx = X, s = model$best_lambda,
                                 type = type, newoffset = newoffset, ...)
  as.vector(yhat)
}
