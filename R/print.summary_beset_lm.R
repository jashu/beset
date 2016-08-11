#' @export

print.summary_beset_lm <- function(object){

  cat(paste("\n=======================================================",
            "\nBest Model:\n",
            object$form,
            "\n\nCoefficients:\n"))
  coefs <- round(object$coef, 3)
  coefs[,4] <- ifelse(coefs[,4] < .001, "<.001", coefs[,4])
  print(coefs, quote = FALSE)
  if(nrow(coefs) > 1){
    cat(paste("\n\nWithin-fold R-squared:", round(object$train_R2,2),
              "\nOut-of-fold R-squared:", round(object$cv_R2,2),
              "\nIndependent R-squared:", round(object$test_R2,2)))
  } else {
    cat("\n\nIntercept-only model: no reliable predictors.")
  }

  cat(paste("\n=======================================================",
            "\nBest Parsimonious Model (within 1 SE of Best Model):\n",
            object$form_1SE,
            "\n\nCoefficients:\n"))
  coefs <- round(object$coef_1SE, 3)
  coefs[,4] <- ifelse(coefs[,4] < .001, "<.001", coefs[,4])
  print(coefs, quote = FALSE)
  if(nrow(coefs) > 1){
    cat(paste("\n\nWithin-fold R-squared:", round(object$train_R2_1SE,2),
              "\nOut-of-fold R-squared:", round(object$cv_R2_1SE,2),
              "\nIndependent R-squared:", round(object$test_R2_1SE,2)))
  } else {
    cat("\n\nIntercept-only model: no reliable predictors.")
  }
  cat("\n=======================================================")
}

#' @export

print.beset_lm <- function(object) summary(object)
