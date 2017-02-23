#' @export

print.summary_beset_elnet <- function(object){
  if(!is.null(object$cv_R2)){

    if(object$best_alpha < .25){
      info <- paste("Primarily ridge with lambda =", object$best_lambda)
    } else if (object$best_alpha > .75){
      info <- paste("Primarily lasso with lambda =", object$best_lambda)
    } else {
      info <- paste("Mixture of ridge and lasso with lambda =",
                    object$best_lambda)
    }
    cat(paste("\n=======================================================",
              "\nBest Model:\n", info,
              "\n\nNon-zero coefficients ranked in order of importance:\n"))
    best_coef <- dplyr::select(object$var_imp, variable, coef = best)
    best_coef <- dplyr::filter(best_coef, abs(coef) > 0)
    best_coef <- dplyr::arrange(best_coef, desc(abs(coef)))
    best_coef <- dplyr::mutate(best_coef, coef = round(coef, 3))
    best_coef <- as.data.frame(best_coef)
    if(nrow(best_coef) > 1){
      print(best_coef, quote = FALSE)
      cat(paste("\n\nTraining  R-squared:", round(object$train_R2,2),
                "\nCross-val R-squared:", round(object$cv_R2,2)))
    } else {
      cat("\n\nNo reliable predictors.")
    }

  }
  if(!is.null(object$cv_R2_1SE)){
    if(object$best_alpha_1SE < .25){
      info <- paste("Primarily ridge with lambda =", object$best_lambda_1SE)
    } else if (object$best_alpha_1SE > .75){
      info <- paste("Primarily lasso with lambda =", object$best_lambda_1SE)
    } else {
      info <- paste("Mixture of ridge and lasso with lambda =",
                    object$best_lambda_1SE)
    }
    cat(paste("\n=======================================================",
              "\nMost Parsimonious Model (within 1 SE of Best Model):\n",
              info,
              "\n\nNon-zero coefficients ranked in order of importance:\n"))
    best_coef <- dplyr::select(object$var_imp, variable, coef = best_sparse)
    best_coef <- dplyr::filter(best_coef, abs(coef) > 0)
    best_coef <- dplyr::arrange(best_coef, desc(abs(coef)))
    best_coef <- dplyr::mutate(best_coef, coef = round(coef, 3))
    best_coef <- as.data.frame(best_coef)
    if(nrow(best_coef) > 1){
      print(best_coef, quote = FALSE)
      cat(paste("\n\nTraining  R-squared:", round(object$train_R2_1SE,2),
                "\nCross-val R-squared:", round(object$cv_R2_1SE,2)))
    } else {
      cat("\n\nNo reliable predictors.")
    }
    cat("\n=======================================================")
  }
}

#' @export
print.beset_elnet <- function(object) summary(object)

#' @export
print.cv_R2 <- function(object)
  cat(paste("Cross-validated R-squared = ", round(object$cv_R2,2), ", 95% CI [",
        round(object$`95% CI`[1],2), ", ", round(object$`95% CI`[2],2), "]",
        sep = ""))

