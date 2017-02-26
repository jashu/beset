#' Summary Method for the \code{beset_glm} Class
#'
#' \code{summary.beset_glm} summarizes the output of a \code{\link{beset_glm}}
#'  object.
#'
#' @param object An object of class \code{\link{beset_glm}}.
#'
#' @export

summary.beset_glm <- function(object, n_pred = NULL, metric = "MCE",
                              oneSE = TRUE, n_cores = 2){
  metric <- match.arg(metric, c("AIC", "MCE", "MSE", "R2"))
  if(is.na(metric)) stop("invalid 'metric' argument")
  best_model <- object$best_AIC
  best_form <- object$fit_stats$form[1]
  if(metric != "AIC"){
    best_metric <- switch(
      metric,
      MCE = min(object$xval_stats$MCE, na.rm = T)[1],
      MSE = min(object$xval_stats$MSE, na.rm = T)[1],
      R2 = max(object$xval_stats$R2, na.rm = T)[1]
    )
    best_idx <- switch(
      metric,
      MCE = which.min(object$xval_stats$MCE),
      MSE = which.min(object$xval_stats$MSE),
      R2 = which.max(object$xval_stats$R2)
    )
    best_form <- object$xval_stats$form[best_idx]
    if(oneSE){
      boundary <- switch(
        metric,
        MCE = best_metric + object$xval_stats$MCE_SE[best_idx],
        MSE = best_metric + object$xval_stats$MSE_SE[best_idx],
        R2 = best_metric - object$xval_stats$R2_SE[best_idx])
      best_1SE <- switch(
        metric,
        MCE = object$xval_stats[object$xval_stats$MCE < boundary,],
        MSE = object$xval_stats[object$xval_stats$MSE < boundary,],
        R2 = object$xval_stats[object$xval_stats$R2 > boundary,]
      )
      best_form <- best_1SE$form[which.min(best_1SE$n_pred)]
    }
    best_model <- if(class(object$best_AIC)[1] == "negbin"){
      update(object$best_AIC, best_form, data = object$model_data)
    } else {
      update(object$best_AIC, best_form,
             family = object$best_AIC$family,
             data = object$model_data)
    }
  }
  best <- c(summary(best_model), loglik = logLik(best_model))
  R2 <- object$fit_stats$R2[object$fit_stats$form == best_form]
  R2_test <- object$test_stats$R2[object$test_stats$form == best_form]
  R2_cv <- do.call("cv_r2",
                   args = c(list(object = best_model, n_cores = n_cores),
                            object$xval_params))

  structure(list(best = best, best_form = best_form, R2 = R2,
                 R2_cv = R2_cv, R2_test = R2_test), class = "summary_beset_glm")
}
