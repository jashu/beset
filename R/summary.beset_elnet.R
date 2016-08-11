#' Summary Method for the \code{beset_lm} Class
#'
#' \code{summary.beset_lm} summarizes the output of a \code{\link{beset_lm}}
#'  object.
#'
#' @param object An object of class \code{\link{beset_lm}}.
#'
#' @import glmnet
#' @export

summary.beset_elnet <- function(object){
  train_R2 <- NULL; cv_R2 <- NULL
  train_R2_1SE <- NULL; cv_R2_1SE <- NULL
  if(!is.null(object$best_model)){
    vars <- rownames(coef(object$best_model, s = object$best_lambda))[-1]
  } else {
    vars <- rownames(coef(object$best_model_1SE,
                          s = object$best_lambda_1SE))[-1]
  }
  var_imp <- dplyr::data_frame(variable = vars, best = 0, best_sparse = 0)
  if(!is.null(object$best_model)){
    coefs <- coef(object$best_model, s = object$best_lambda)
    i <- coefs@i[-1]
    var_imp$best[i] <- coefs@x[-1]
    best_lambda <- object$best_model$lambda
    if(min(best_lambda) > object$best_lambda)
      object$best_lambda <- min(best_lambda)
    if(max(best_lambda) < object$best_lambda)
      object$best_lambda <- max(best_lambda)
    train_R2 <- approx(object$best_model$lambda, object$best_model$dev.ratio,
                       object$best_lambda)$y
    results <- object$results[object$results$alpha == object$best_alpha,]
    null_dev <- object$best_model$nulldev
    cv_dev <- approx(results$lambda, results$cve,
                     object$best_lambda)$y * object$best_model$nobs
    cv_R2 <- 1 - cv_dev / null_dev
  }
  if(!is.null(object$best_model_1SE)){
    coefs <- coef(object$best_model_1SE, s = object$best_lambda_1SE)
    i <- coefs@i[-1]
    var_imp$best_sparse[i] <- coefs@x[-1]
    best_lambda <- object$best_model_1SE$lambda
    if(min(best_lambda) > object$best_lambda_1SE)
      object$best_lambda_1SE <- min(best_lambda)
    if(max(best_lambda) < object$best_lambda_1SE)
      object$best_lambda_1SE <- max(best_lambda)
    train_R2_1SE <- approx(object$best_model_1SE$lambda,
                           object$best_model_1SE$dev.ratio,
                           object$best_lambda_1SE)$y
    results <- object$results[object$results$alpha == object$best_alpha_1SE,]
    null_dev <- object$best_model_1SE$nulldev
    cv_dev <- approx(results$lambda, results$cve,
                     object$best_lambda_1SE)$y * object$best_model_1SE$nobs
    cv_R2_1SE <- 1 - cv_dev / null_dev
  }
  var_imp <- dplyr::arrange(var_imp, desc(abs(best + best_sparse)))
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
