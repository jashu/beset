#' @export
print.R2 <- function(object){
  cat(paste("Model-fit R-squared =", round(object$R2fit,2)))
  if(!is.null(object$R2new)){
    cat(paste(", Predictive R-squared =", round(object$R2new, 2)))
  }
  if(!is.null(object$R2cv)){
    cat("\n")
    print(object$R2cv)
  }
}
#' @export
print.cv_R2 <- function(object)
  cat(paste("Cross-validated R-squared = ", round(object$cv_R2,2), ", 95% CI [",
            round(object$`95% CI`[1],2), ", ", round(object$`95% CI`[2],2), "]",
            sep = ""))
#' @export
print.beset_glm <- function(object) print(summary(object))
#' @export
print.beset_zeroinfl <- function(object) print(summary(object))
#' @export
print.summary_beset_glm <- function(
  object, digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"), ...){
  cat("\n=======================================================",
      "\nBest Model:\n ", object$best_form, "\n")
  if(length(object$near_best) > 0){
    cat("\nNearly Equivalent Models:", object$near_best, sep = "\n  ")
  }
  x <- object$best
  cat("\nDeviance Residuals: \n")
  if (x$df.residual > 5) {
    x$deviance.resid <- setNames(quantile(x$deviance.resid,
                                          na.rm = TRUE),
                                 c("Min", "1Q", "Median", "3Q", "Max"))
  }
  xx <- zapsmall(x$deviance.resid, digits + 1L)
  print.default(xx, digits = digits, na.print = "", print.gap = 2L)
  if (length(x$aliased) == 0L) {
    cat("\nNo Coefficients\n")
  }
  else {
    df <- if ("df" %in% names(x))
      x[["df"]]
    else NULL
    if (!is.null(df) && (nsingular <- df[3L] - df[1L]))
      cat("\nCoefficients: (", nsingular,
          " not defined because of singularities)\n",
          sep = "")
    else cat("\nCoefficients:\n")
    coefs <- x$coefficients
    if (!is.null(aliased <- x$aliased) && any(aliased)) {
      cn <- names(aliased)
      coefs <- matrix(NA, length(aliased), 4L, dimnames = list(cn,
                                                               colnames(coefs)))
      coefs[!aliased, ] <- x$coefficients
    }
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
                 na.print = "NA", ...)
  }
  cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ",
      format(x$dispersion), ")\n\n", sep = "")
  if (nzchar(mess <- naprint(x$na.action)))
    cat("  (", mess, ")\n", sep = "")
  cat("Log-likelihood: ", formatC(x$loglik, digits = digits), " on ",
      x$loglik_df, " Df\nAIC: ",
      format(x$aic, digits = max(4L, digits + 1L)),
      "\n\n", "Number of Fisher Scoring iterations: ", x$iter,
      "\n", sep = "")
  correl <- x$correlation
  if (!is.null(correl)) {
    p <- NCOL(correl)
    if (p > 1) {
      cat("\nCorrelation of Coefficients:\n")
      correl <- format(round(correl, 2L), nsmall = 2L,
                       digits = digits)
      correl[!lower.tri(correl)] <- ""
      print(correl[-1, -p, drop = FALSE], quote = FALSE)
    }
  }
  if("summary.negbin" %in% class(x)){
    dp <- max(2 - floor(log10(x$SE.theta)), 0)
    cat("\n              Theta: ", format(round(x$theta, dp), nsmall = dp),
        "\n          Std. Err.: ", format(round(x$SE.theta, dp), nsmall = dp),
        "\n")
    if (!is.null(x$th.warn))
      cat("Warning while fitting theta:", x$th.warn, "\n")
  }
  cat("\n")
  cat(paste("Train-sample R-squared =", round(object$R2,2)))
  if(!is.null(object$R2_test)){
    cat(paste(", Test-sample R-squared =", round(object$R2_test, 2)))
  }
  cat("\n")
  print(object$R2_cv)
  cat("\n=======================================================")
}
#' @export
print.summary_beset_zeroinfl <- function(
  object, digits = max(3, getOption("digits") - 3), ...){
  cat(paste("\n=======================================================",
            "\nBest Model:\n",
            object$best_form,
            "\n\n"))
  x <- object$best
  if (!x$converged) {
    cat("model did not converge\n")
  } else {
    cat("Pearson residuals:\n")
    print(structure(quantile(x$residuals),
                    names = c("Min","1Q", "Median", "3Q", "Max")),
          digits = digits, ...)
    cat(paste("\nCount model coefficients (", x$dist, " with log link):\n",
              sep = ""))
    printCoefmat(x$coefficients$count, digits = digits, signif.legend = FALSE)
    cat(paste("\nZero-inflation model coefficients (binomial with ",
              x$link, " link):\n", sep = ""))
    printCoefmat(x$coefficients$zero, digits = digits, signif.legend = FALSE)
    if (getOption("show.signif.stars") &
        any(rbind(x$coefficients$count, x$coefficients$zero)[, 4] < 0.1))
      cat("---\nSignif. codes: ",
          "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1",
          "\n")
    if (x$dist == "negbin")
      cat(paste("\nTheta =", round(x$theta, digits), "\n")) else cat("\n")
    cat(paste("Number of iterations in", x$method, "optimization:",
              tail(na.omit(x$optim$count), 1), "\n"))
    cat("Log-likelihood:", formatC(x$loglik, digits = digits),
        "on", x$n - x$df.residual, "Df\nAIC: ",
        format(x$aic, digits = max(4L, digits + 1L)), "\n\n")
    cat(paste("Train-sample R-squared =", round(object$R2,2)))
    if(!is.null(object$R2test)){
      cat(paste(", Test-sample R-squared =", round(object$R2_test, 2)))
    }
    cat("\n")
    print(object$R2_cv)
  }
  cat("\n=======================================================")
}

#' @export
print.summary_beset_elnet <- function(object){
  cat("\n=======================================================",
      "\nBest Model:\n", sep = "")
  if(object$best_alpha < .25){
    cat("Primarily ridge ")
  } else if (object$best_alpha > .75){
    cat("Primarily lasso ")
  } else {
    cat("Mixture of ridge and lasso ")
  }
  cat("(alpha = ", object$best_alpha, ") with lambda = ", object$best_lambda,
      sep = "")
  cat("\n\nNon-zero coefficients ranked in order of importance:\n")
  best_coef <- dplyr::select(object$var_imp, variable, stand.coef = coef)
  best_coef <- dplyr::filter(best_coef, abs(stand.coef) > 0)
  best_coef <- dplyr::arrange(best_coef, desc(abs(stand.coef)))
  best_coef <- dplyr::mutate(best_coef, stand.coef = round(stand.coef, 3))
  best_coef <- as.data.frame(best_coef)
  if(nrow(best_coef) > 1){
    print(best_coef, quote = FALSE)
    cat(paste("\nTrain-sample R-squared =", round(object$R2,2)))
    if(!is.null(object$R2test)){
      cat(paste(", Test-sample R-squared =", round(object$R2_test, 2)))
    }
    cat("\n")
    print(object$R2_cv)
  } else {
    cat("\n\nNo reliable predictors.")
  }
  cat("\n=======================================================")
}

#' @export
print.beset_elnet <- function(object) print(summary(object))





