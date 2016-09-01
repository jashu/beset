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
  best <- summary(object$best_model)
  form <- as.character(terms(best))
  form <- paste0(form[2], " ", form[1], " ", form[3])
  best_form <- form
  best_cve <- data$cv_CE[data$form == best_form]
  best_1SE <- summary(object$best_model_1SE)
  form <- as.character(terms(best_1SE))
  form <- paste0(form[2], " ", form[1], " ", form[3])
  best_form_1SE <- form
  best_cve_1SE <- data$cv_CE[data$form == best_form_1SE]
  structure(list(best = best, best_1SE = best_1SE,
                 best_form = best_form, best_form_1SE = best_form_1SE,
                 best_cve = best_cve, best_cve_1SE = best_cve_1SE),
            class = "summary_beset_glm")
}
