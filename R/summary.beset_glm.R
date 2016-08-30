#' Summary Method for the \code{beset_glm} Class
#'
#' \code{summary.beset_glm} summarizes the output of a \code{\link{beset_glm}}
#'  object.
#'
#' @param object An object of class \code{\link{beset_glm}}.
#'
#' @export

summary.beset_glm <- function(object){
  data <- object$best_subsets
  n_pred <- length(all.vars(object$best_model$terms)) - 1
  form <- data$form[data$n_pred == n_pred]
  coef <- summary(object$best_model)$coefficients
  train_R2 <- data$train_R2[data$n_pred == n_pred]
  cv_R2 <- data$cv_R2[data$n_pred == n_pred]
  test_R2 <- data$test_R2[data$n_pred == n_pred]
  n_pred <- length(all.vars(object$best_model_1SE$terms)) - 1
  form_1SE <- data$form[data$n_pred == n_pred]
  coef_1SE <- summary(object$best_model_1SE)$coefficients
  train_R2_1SE<- data$train_R2[data$n_pred == n_pred]
  cv_R2_1SE <- data$cv_R2[data$n_pred == n_pred]
  test_R2_1SE <- data$test_R2[data$n_pred == n_pred]
  structure(list(form = form,
                 coef = coef,
                 train_R2 = train_R2,
                 cv_R2 = cv_R2,
                 test_R2 = test_R2,
                 form_1SE = form_1SE,
                 coef_1SE = coef_1SE,
                 train_R2_1SE = train_R2_1SE,
                 cv_R2_1SE = cv_R2_1SE,
                 test_R2_1SE = test_R2_1SE),
            class = "summary_beset_glm")
}
