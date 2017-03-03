#' Summary Method for the \code{beset_glm} Class
#'
#' \code{summary.beset_glm} summarizes the output of a \code{\link{beset_glm}}
#'  object.
#'
#' @param object An object of class \code{\link{beset_glm}}.
#'
#' @export

summary.beset_glm <- function(object, metric = "MCE", n_pred = NULL,
                              oneSE = TRUE, n_cores = 2){
  metric <- tryCatch(match.arg(metric, c("AIC", "MCE", "MSE", "R2")),
                     error = function(c){
                       c$message <- gsub("arg", "metric", c$message)
                       c$call <- NULL
                       stop(c)
                     })
  best_model <- object$best_AIC
  best_form <- object$stats$fit$form[1]
  if(!is.null(n_pred)){
    best_form <- object$stats$cv$form[object$stats$cv$n_pred == n_pred]
    best_model <- if(class(best_model)[1] == "negbin"){
      update(best_model, best_form, data = object$model_data)
    } else {
      update(best_model, best_form, family = best_model$family,
             data = object$model_data)
    }
  }
  if(metric != "AIC" && is.null(n_pred)){
    best_metric <- switch(
      metric,
      MCE = min(object$stats$cv$MCE, na.rm = T)[1],
      MSE = min(object$stats$cv$MSE, na.rm = T)[1],
      R2 = max(object$stats$cv$R2, na.rm = T)[1]
    )
    best_idx <- switch(
      metric,
      MCE = which.min(object$stats$cv$MCE),
      MSE = which.min(object$stats$cv$MSE),
      R2 = which.max(object$stats$cv$R2)
    )
    best_form <- object$stats$cv$form[best_idx]
    if(oneSE){
      boundary <- switch(
        metric,
        MCE = best_metric + object$stats$cv$MCE_SE[best_idx],
        MSE = best_metric + object$stats$cv$MSE_SE[best_idx],
        R2 = best_metric - object$stats$cv$R2_SE[best_idx])
      best_1SE <- switch(
        metric,
        MCE = object$stats$cv[object$stats$cv$MCE < boundary,],
        MSE = object$stats$cv[object$stats$cv$MSE < boundary,],
        R2 = object$stats$cv[object$stats$cv$R2 > boundary,]
      )
      best_form <- best_1SE$form[which.min(best_1SE$n_pred)]
    }
    best_model <- if(class(best_model)[1] == "negbin"){
      update(best_model, best_form, data = object$model_data)
    } else {
      update(best_model, best_form, family = best_model$family,
             data = object$model_data)
    }
  }
  p <- object$stats$fit$n_pred[object$stats$fit$form == best_form]
  loglik <- logLik(best_model)
  best <- c(summary(best_model), loglik = loglik,
            loglik_df = attr(loglik, "df"))
  R2 <- object$stats$fit$R2[object$stats$fit$form == best_form]
  R2_test <- object$stats$test$R2[object$stats$test$form == best_form]
  R2_cv <- do.call("cv_r2",
                   args = c(list(object = best_model, n_cores = n_cores),
                            object$cv_params))
  near_equals <- dplyr::filter(object$stats$fit, n_pred == p)
  if(metric == "AIC") metric <- "MCE"
  near_equals <- dplyr::select(near_equals, form,
                               metric = dplyr::starts_with(metric))
  if(metric == "MCE") near_equals$metric <- exp(-near_equals$metric)
  best_metric <- if(metric == "MSE") min(near_equals$metric, na.rm = T)[1] else
    max(near_equals$metric, na.rm = T)[1]
  boundary <- c(best_metric - .01*best_metric, best_metric + .01*best_metric)
  near_best <- near_equals$form[dplyr::between(
    near_equals$metric, boundary[1], boundary[2])]
  near_best <- near_best[near_best != best_form]

  structure(list(best = best, best_form = best_form, near_best = near_best,
                 R2 = R2, R2_cv = R2_cv, R2_test = R2_test),
            class = "summary_beset_glm")
}
