#' Predicted R-squared
#'
#' \code{predict_R2} calculates predicted R-squared, which indicates how well a
#' regression model predicts responses for new observations.
#'
#' \code{predict_R2} provides a metric analagous to R-squared, except instead of
#' describing how well a linear model fits the original data, it describes how
#' well the model predicts new data. If new data are not provided then the
#' predicted R-square will be estimated using the PRESS statistic
#' (predicted residual sum of squares), which is given by the sum of squares of
#' the ratio of the model residuals to the leverages (the diagonal elements of
#' the hat matrix) subtracted from 1. This provides for fast computation of the
#' leave-one-out residuals from the data returned by a single call to
#' \code{\link[base]{lm}}, without needing to recompute the model N times after
#' removing each of N observations. Note that predicted R-squared can be negative
#' if overfitting is severe; i.e., when predicting new data, the model performs
#' worse than if one were to always predict the sample mean.
#'
#' @param object Object of class inheriting from "lm"
#'
#' @param newdata	An optional data frame in which to look for variables with
#' which to predict. If omitted, leave-one-out cross validation is used to
#' simulate predictions on new data.
#'
#' @return Predicted R-squared as 1 - RSS / TSS, where RSS corresponds to the
#' sum of the squared residuals (observed y - predicted y) and TSS is the total
#' sum of squares.
#'
#' @export

predict_R2 <- function(object, newdata = NULL){
  if(!inherits(object, "lm"))
    stop("Object must be a linear model.")
  if(is.null(newdata)){
    pRSS <- sum((residuals(object) / (1 - influence(object)$hat))^2)
    TSS <- sum(anova(object)[2])
  }else{
    y_obs <- unlist(newdata[, names(object$model)[1]])
    y_bar <- mean(y_obs, na.rm = TRUE)
    y_hat <- predict(object, newdata)
    pRSS <- sum((y_obs - y_hat)^2, na.rm = TRUE)
    TSS <- sum((y_obs - y_bar)^2, na.rm = TRUE)
  }
  1 - pRSS/TSS
}
