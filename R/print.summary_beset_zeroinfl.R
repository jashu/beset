#' @export

print.summary_beset_zeroinfl <- function(object,
                                         digits = max(3,
                                                      getOption("digits") - 3),
                                         ...){
  cat(paste("\n=======================================================",
            "\nBest Model:\n",
            object$best_form,
            "\n\n"))
  x <- object$best
  if (!x$converged) {
    cat("model did not converge\n")
  }
  else {
    cat(paste("Count model coefficients (", x$dist, " with log link):\n",
              sep = ""))
    print.default(format(x$coefficients$count, digits = digits),
                  print.gap = 2, quote = FALSE)
    if (x$dist == "negbin")
      cat(paste("Theta =", round(x$theta, digits), "\n"))
    cat(paste("\nZero-inflation model coefficients (binomial with ",
              x$link, " link):\n", sep = ""))
    print.default(format(x$coefficients$zero, digits = digits),
                  print.gap = 2, quote = FALSE)
  }
  cat(paste("\n\nTraining R-squared:", round(object$train_R2,2),
              "\nCross-validation R-squared:", round(object$cv_R2,2),
              "\nTest R-squared:", round(object$test_R2,2)))
  cat(paste("\n=======================================================",
            "\nBest Parsimonious Model (within 1 SE of Best Model):\n",
            object$best_form_1SE,
            "\n\n"))
  x <- object$best_1SE
  if (!x$converged) {
    cat("model did not converge\n")
  }
  else {
    cat(paste("Count model coefficients (", x$dist, " with log link):\n",
              sep = ""))
    print.default(format(x$coefficients$count, digits = digits),
                  print.gap = 2, quote = FALSE)
    if (x$dist == "negbin")
      cat(paste("Theta =", round(x$theta, digits), "\n"))
    cat(paste("\nZero-inflation model coefficients (binomial with ",
              x$link, " link):\n", sep = ""))
    print.default(format(x$coefficients$zero, digits = digits),
                  print.gap = 2, quote = FALSE)
  }
    cat(paste("\n\nTraining R-squared:", round(object$train_R2_1SE,2),
              "\nCross-validation R-squared:", round(object$cv_R2_1SE,2),
              "\nTest R-squared:", round(object$test_R2_1SE,2)))
  cat("\n=======================================================")
}

#' @export
print.beset_zeroinfl <- function(object) print(summary(object))
