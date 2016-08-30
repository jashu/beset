#' Summary Method for the \code{beset_zeroinfl} Class
#'
#' \code{summary.beset_zeroinfl} summarizes the output of a
#' \code{\link{beset_zeroinfl}} object.
#'
#' @param object An object of class \code{\link{beset_zeroinfl}}.
#'
#' @export

summary.beset_zeroinfl <- function(object){
  data <- object$best_subsets
  best <- summary(object$best_model)
  form <- as.character(best$formula)
  form <- paste0(form[2], " ", form[1], " ", form[3])
  best_form <- form
  idx <- which(data$form == form)
  train_R2 <- data$train_R2[idx]
  cv_R2 <- data$cv_R2[idx]
  test_R2 <- data$test_R2[idx]
  best_1SE <- summary(object$best_model_1SE)
  form <- as.character(best_1SE$formula)
  form <- paste0(form[2], " ", form[1], " ", form[3])
  best_form_1SE <- form
  idx <- which(data$form == form)
  train_R2_1SE <- data$train_R2[idx]
  cv_R2_1SE <- data$cv_R2[idx]
  test_R2_1SE <- data$test_R2[idx]

  structure(list(best = best, best_1SE = best_1SE,
                 best_form = best_form, best_form_1SE = best_form_1SE,
                 train_R2 = train_R2, train_R2_1SE = train_R2_1SE,
                 cv_R2 = cv_R2, cv_R2_1SE = cv_R2_1SE,
                 test_R2 = test_R2, test_R2_1SE = test_R2_1SE),
            class = "summary_beset_zeroinfl")
}

