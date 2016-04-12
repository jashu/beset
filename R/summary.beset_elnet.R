#' Summary Method for the \code{beset_lm} Class
#'
#' \code{summary.beset_lm} summarizes the output of a \code{\link{beset_lm}}
#'  object.
#'
#' @param object An object of class \code{\link{beset_lm}}.
#'
#' @export

summary.beset_elnet <- function(object){
  best_coefs <- coef(object$best_model, s = object$best_lambda)
  i <- best_coefs@i[-1]
  vars <- best_coefs@Dimnames[[1]][-1]
  var_imp <- dplyr::data_frame(variable = vars, best = 0, best_sparse = 0)
  var_imp$best[i] <- best_coefs@x[-1]
  best_coefs <- coef(object$best_model_1SE, s = object$best_lambda_1SE)
  i <- best_coefs@i[-1]
  var_imp$best_sparse[i] <- best_coefs@x[-1]
  var_imp <- dplyr::arrange(var_imp, desc(abs(best + best_sparse)))
  train_R2 <- approx(object$best_model$lambda, object$best_model$dev.ratio,
                     object$best_lambda)$y
  train_R2_1SE <- approx(object$best_model_1SE$lambda,
                         object$best_model_1SE$dev.ratio,
                         object$best_lambda_1SE)$y
  results <- object$results[object$results$alpha == object$best_alpha,]
  null_dev <- object$best_model$nulldev
  cv_dev <- approx(results$lambda, results$cve,
                   object$best_lambda)$y * object$best_model$nobs
  cv_R2 <- 1 - cv_dev / null_dev
  results <- object$results[object$results$alpha == object$best_alpha_1SE,]
  null_dev <- object$best_model_1SE$nulldev
  cv_dev <- approx(results$lambda, results$cve,
                   object$best_lambda_1SE)$y * object$best_model_1SE$nobs
  cv_R2_1SE <- 1 - cv_dev / null_dev

  ans <- structure(list(var_imp = var_imp,
                        train_R2 = train_R2,
                        cv_R2 = cv_R2,
                        train_R2_1SE = train_R2_1SE,
                        cv_R2_1SE = cv_R2_1SE,
                        best_alpha = object$best_alpha,
                        best_alpha_1SE = object$best_alpha_1SE,
                        best_lambda = object$best_lambda,
                        best_lambda_1SE = object$best_lambda_1SE),
                   class = "summary_beset_elnet")
  return(ans)
}
