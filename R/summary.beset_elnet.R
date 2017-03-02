#' Summary Method for the \code{beset_elnet} Class
#'
#' \code{summary.beset_elnet} summarizes the output of a
#' \code{\link{beset_elnet}} object.
#'
#' @param object An object of class \code{\link{beset_elnet}}.
#'
#' @import glmnet
#' @export

summary.beset_elnet <- function(object, metric = "MCE", oneSE = TRUE,
                                n_cores = 2){
  metric <- match.arg(metric, c("MCE", "MSE", "R2"))
  if(is.na(metric)) stop("invalid 'metric' argument")
  best_metric <- switch(
      metric,
      MCE = min(object$stats$xval$MCE, na.rm = T)[1],
      MSE = min(object$stats$xval$MSE, na.rm = T)[1],
      R2 = max(object$stats$xval$R2, na.rm = T)[1]
    )
    best_idx <- switch(
      metric,
      MCE = which.min(object$stats$xval$MCE),
      MSE = which.min(object$stats$xval$MSE),
      R2 = which.max(object$stats$xval$R2)
    )
  best_alpha <- object$stats$xval$alpha[best_idx]
  best_lambda <- object$stats$xval$lambda[best_idx]
  if(oneSE){
    boundary <- switch(
      metric,
      MCE = best_metric + object$stats$xval$MCE_SE[best_idx],
      MSE = best_metric + object$stats$xval$MSE_SE[best_idx],
      R2 = best_metric - object$stats$xval$R2_SE[best_idx])
    best_1SE <- switch(
      metric,
      MCE = object$stats$xval[object$stats$xval$MCE < boundary,],
      MSE = object$stats$xval[object$stats$xval$MSE < boundary,],
      R2 = object$stats$xval[object$stats$xval$R2 > boundary,]
    )
    best_alpha <- max(best_1SE$alpha)
    best_1SE <- dplyr::filter(best_1SE, alpha == best_alpha)
    best_lambda <- max(best_1SE$lambda)
    best_idx <- which(object$stats$fit$alpha == best_alpha &
                        object$stats$fit$lambda == best_lambda)
  }
  best_model <- object$model_fits[[as.character(best_alpha)]]
  coefs <- coef(best_model, s = best_lambda)
  var_imp <- dplyr::data_frame(variable = rownames(coefs)[-1], coef = 0)
  var_imp$coef[coefs@i[-1]] <- coefs@x[-1]
  var_imp <- dplyr::arrange(var_imp, dplyr::desc(abs(coef)))
  R2 <- object$stats$fit$R2[best_idx]
  R2_test <- object$stats$test$R2[best_idx]
  x <- object$model_params$x
  y <- object$model_params$y
  family <- object$model_params$family
  cv_R2 <- sapply(object$xval_params$fold_ids, function(i){
    fit <- update(best_model, x = x[i,], y = y[i], family = family,
                  alpha = best_alpha)
    y_hat <- predict(fit, x[-i,], best_lambda, "response")
    predict_metrics_(y[-i], y_hat, family)$R_squared
  })
  R2_boot <- boot::boot(cv_R2, function(x,i) median(x[i]), 1000)
  boot_CI <- boot::boot.ci(R2_boot, type = "bca")
  R2_cv <- structure(
    list(cv_R2 = median(cv_R2, na.rm = TRUE), `95% CI` = boot_CI$bca[4:5],
                 R2 = cv_R2), class = "cv_R2")

  structure(list(best = best_model, best_alpha = best_alpha,
                 best_lambda = best_lambda, var_imp = var_imp,
                 R2 = R2, R2_cv = R2_cv, R2_test = R2_test),
            class = "summary_beset_elnet")
}
