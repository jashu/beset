#' @import gcdnet
#' @export

predict.rr_elnet <- function(object, data){
  data <- data[, rownames(object$cv_gcdnet$gcdnet.fit$beta)]
  predict(object$cv_gcdnet, data)
}
