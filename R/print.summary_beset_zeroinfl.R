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
  } else {
    cat("Pearson residuals:\n")
    print(structure(quantile(x$residuals), names = c("Min",
                                                     "1Q", "Median", "3Q", "Max")), digits = digits, ...)
    cat(paste("\nCount model coefficients (", x$dist, " with log link):\n",
              sep = ""))
    printCoefmat(x$coefficients$count, digits = digits, signif.legend = FALSE)
    cat(paste("\nZero-inflation model coefficients (binomial with ",
              x$link, " link):\n", sep = ""))
    printCoefmat(x$coefficients$zero, digits = digits, signif.legend = FALSE)
    if (getOption("show.signif.stars") & any(rbind(x$coefficients$count,
                                                   x$coefficients$zero)[, 4] < 0.1))
      cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1",
          "\n")
    if (x$dist == "negbin")
      cat(paste("\nTheta =", round(x$theta, digits), "\n")) else cat("\n")
    cat(paste("Number of iterations in", x$method, "optimization:",
              tail(na.omit(x$optim$count), 1), "\n"))
    cat("Log-likelihood:", formatC(x$loglik, digits = digits),
        "on", x$n - x$df.residual, "Df\n")
    cat("Cross-validation error:", round(object$best_cve, 2),
        "(cross entropy)\n")
  }
  cat(paste("\n=======================================================",
            "\nBest Parsimonious Model (within 1 SE of Best Model):\n",
            object$best_form_1SE,
            "\n\n"))
  if(identical(object$best, object$best_1SE)){
    cat("identical to Best Model (see above)")
  } else {
    x <- object$best_1SE
    if (!x$converged) {
      cat("model did not converge\n")
    }
    else {
      cat("Pearson residuals:\n")
      print(structure(quantile(x$residuals), names = c("Min",
                                                       "1Q", "Median", "3Q", "Max")), digits = digits, ...)
      cat(paste("\nCount model coefficients (", x$dist, " with log link):\n",
                sep = ""))
      printCoefmat(x$coefficients$count, digits = digits, signif.legend = FALSE)
      cat(paste("\nZero-inflation model coefficients (binomial with ",
                x$link, " link):\n", sep = ""))
      printCoefmat(x$coefficients$zero, digits = digits, signif.legend = FALSE)
      if (getOption("show.signif.stars") & any(rbind(x$coefficients$count,
                                                     x$coefficients$zero)[, 4] < 0.1))
        cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1",
            "\n")
      if (x$dist == "negbin")
        cat(paste("\nTheta =", round(x$theta, digits), "\n")) else cat("\n")
      cat(paste("Number of iterations in", x$method, "optimization:",
                tail(na.omit(x$optim$count), 1), "\n"))
      cat("Log-likelihood:", formatC(x$loglik, digits = digits),
          "on", x$n - x$df.residual, "Df\n")
      cat("Cross-validation error:", round(object$best_cve_1SE, 2),
          "(cross entropy)\n")
    }
  }
  cat("\n=======================================================")
}

#' @export
print.beset_zeroinfl <- function(object) print(summary(object))
