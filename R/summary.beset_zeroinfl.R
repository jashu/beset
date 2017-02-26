#' Summary Method for the \code{beset_zeroinfl} Class
#'
#' \code{summary.beset_zeroinfl} summarizes the output of a
#' \code{\link{beset_zeroinfl}} object.
#'
#' @param object An object of class \code{\link{beset_zeroinfl}}.
#'
#' @export
summary.beset_zeroinfl <- function(object, n_pred = NULL, metric = "MCE",
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
      n_pred <- best_1SE$n_count_pred + best_1SE$n_zero_pred
      best_form <- best_1SE$form[which.min(n_pred)]
    }
    best_model <- pscl::zeroinfl(formula = formula(best_form),
                                 data = object$model_data,
                                 weights = object$best_AIC$weights,
                                 dist = object$best_AIC$dist,
                                 link = object$best_AIC$link,
                                 control = object$best_AIC$control)
  }
  best <- c(summary(best_model), aic = AIC(best_model))
  R2 <- object$fit_stats$R2[object$fit_stats$form == best_form]
  R2_test <- object$test_stats$R2[object$test_stats$form == best_form]
  R2_cv <- do.call("cv_r2",
                   args = c(list(object = best_model, n_cores = n_cores),
                            object$xval_params))

  structure(list(best = best, best_form = best_form, R2 = R2,
                 R2_cv = R2_cv, R2_test = R2_test),
            class = "summary_beset_zeroinfl")
}
