#' Summarizing Bootstrap Objects
#'
#' Methods for the \code{bootstrap} class.
#'
#' @param object An object of class \code{bootstrap}.
#' @export

summary.bootstrap <- function (object, ...){
  summary(object$stats, ...)
}

#' @export
print.bootstrap <- function(object){
  print(summary(object$stats))
}
