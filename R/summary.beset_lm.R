#' Summary Method for the \code{beset_lm} Class
#'
#' \code{summary.beset_lm} summarizes the output of a \code{\link{beset_lm}}
#'  object.
#'
#' @param object An object of class \code{\link{beset_lm}}.
#'
#' @export

summary.beset_lm <- function(object){
  all_vars <- all.vars(object$best_model$terms)
  x <- "1"
  if(length(all_vars) > 1){
    x <- paste0(all_vars[2:length(all_vars)], collapse = " + ")
  }
  y <- all_vars[1]
  form <- paste(y, "~", x)
  coef <- summary(object$best_model)$coefficients
  train_R2 <-
    mean(object$R2$R2_train[object$R2$n_preds == (length(all_vars) - 1)])
  cv_R2 <-
    mean(object$R2$R2_cv[object$R2$n_preds == (length(all_vars) - 1)])
  test_R2 <-
    mean(object$R2$R2_test[object$R2$n_preds == (length(all_vars) - 1)])
  all_vars <- all.vars(object$best_model_1SE$terms)
  x <- "1"
  if(length(all_vars) > 1){
    x <- paste0(all_vars[2:length(all_vars)], collapse = " + ")
  }
  y <- all_vars[1]
  form_1SE <- paste(y, x, sep = " ~ ")
  coef_1SE <- summary(object$best_model_1SE)$coefficients
  train_R2_1SE <-
    mean(object$R2$R2_train[object$R2$n_preds == (length(all_vars) - 1)])
  cv_R2_1SE <-
    mean(object$R2$R2_cv[object$R2$n_preds == (length(all_vars) - 1)])
  test_R2_1SE <-
    mean(object$R2$R2_test[object$R2$n_preds == (length(all_vars) - 1)])
  ans <- structure(list(form = form,
                        coef = coef,
                        train_R2 = train_R2,
                        cv_R2 = cv_R2,
                        test_R2 = test_R2,
                        form_1SE = form_1SE,
                        coef_1SE = coef_1SE,
                        train_R2_1SE = train_R2_1SE,
                        cv_R2_1SE = cv_R2_1SE,
                        test_R2_1SE = test_R2_1SE),
                   class = "summary_beset_lm")
  return(ans)
}
