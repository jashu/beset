#' Predicted R-squared
#'
#' \code{predict_R2} calculates predicted R-squared, which indicates how well a
#' regression model predicts responses for new observations.
#'
#' \code{predict_R2} provides a metric analagous to R-squared, except instead of
#' describing how well a model fits the original training data, it describes how
#' well the model predicts independent test data. Optionally, if the model was
#' fit using ordinary least squares (i.e., with \code{\link[base]{lm}}, the
#' \code{test_data} argument can be omitted, in which case the predicted
#' R-square will be estimated using the PRESS statistic (predicted residual sum
#' of squares), which is given by the sum of squares of the ratio of the model
#' residuals to the leverages (the diagonal elements of the hat matrix)
#' subtracted from 1. This provides for fast computation of the leave-one-out
#' residuals from information already present in the \code{\link[base]{lm}}
#' object, without needing to recompute the model N times after removing each of
#' N observations.
#'
#' Note that predicted R-squared can be negative if overfitting is severe; i.e.,
#' when predicting new data, the model performs worse than if one were to always
#' predict the sample mean.
#'
#' @param object Object of model class "lm", "glm", or "zeroinfl".
#'
#' @param test_data	Data frame containing independent test data. Variable names
#' must be identical to  which to look for variables with which to predict.
#' Optional for \code{\link[base]{lm}} objects. Required for all other model
#' types.
#'
#' @return Predicted R-squared. Except for logistic models, this is calculated
#' as 1 - RSS / TSS, where RSS corresponds to the sum of the squared residuals
#' (observed y - predicted y) and TSS is the total sum of squares. For logistic
#' models, this is calculated as 1 - L_MOD / L_NULL, where L_MOD is the
#' log-likelihood for the model prediction and L_NULL is the log-likelihood for
#' the null model (that every individual belongs to the majority class with
#' probability equal to the observed proportion of individuals in the majority).
#'
#' @export

predict_R2 <- function(object, test_data = NULL){
  if(is.null(test_data)){
    if(any(class(object) != "lm")){
      stop("'test_data' must be included for this model class")
    }
    error_with_model <- sum((residuals(object) / (1 - influence(object)$hat))^2)
    error_without_model <- sum(anova(object)[2])
  } else {
    y_obs <- unlist(test_data[, names(object$model)[1]])
    if(is.factor(y_obs)) y_obs <- as.integer(y_obs) - 1
    y_bar <- mean(y_obs, na.rm = TRUE)
    y_hat <- predict(object, test_data, type = "response")
    if("glm" %in% class(object) && object$family$family == "binomial"){
      error_with_model <- y_obs * log(y_hat) + (1 - y_obs) * log(1 - y_hat)
      error_with_model <- mean(error_with_model, na.rm = T)
      error_without_model <- y_obs * log(y_bar) + (1 - y_obs) * log(1 - y_bar)
      error_without_model <- mean(error_without_model, na.rm = T)
    } else {
      error_with_model <- sum((y_obs - y_hat)^2, na.rm = TRUE)
      error_without_model <- sum((y_obs - y_bar)^2, na.rm = TRUE)
    }
  }
  1 - error_with_model/error_without_model
}
