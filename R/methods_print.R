#' @importFrom utils tail
#' @import purrr
#' @import dplyr

#' @export
print.beset <- function(x, ...) print(summary(x, ...))

#' @export
print.cross_valid <- function(x, ...){
  cat("Mean predictive performance under ")
  if(x$parameters$n_obs == x$parameters$n_folds){
    cat("leave-one-out cross-validation:")
    } else {
      cat(paste(x$parameters$n_folds, "-fold cross-validation\n", sep = ""))
    }
  if(x$parameters$n_reps > 1){
    cat(paste("(with min-max range over ", x$parameters$n_reps, " repetitions)",
            sep = ""))
  }
  cat("\n")
  results_frame <- results_frame <- data_frame(
    Mean =  map_dbl(x$stats, "mean"),
    S.E. = map_dbl(x$stats, "btwn_fold_se")
  )
  if(x$parameters$n_reps > 1){
    results_frame$Min <- map_dbl(x$stats, ~ .x$btwn_rep_range[1])
    results_frame$Max <- map_dbl(x$stats, ~ .x$btwn_rep_range[2])
  }
  results_frame <- dplyr::mutate_all(results_frame, ~ signif(., 3))
  results_frame <- as.data.frame(results_frame)
  metrics <- names(x$stats)
  if(x$parameters$family != "gaussian") metrics[4] <- "r2d"
  row.names(results_frame) <- map(
    metrics, ~ switch(.x,
                      rsq = "Variance Explained",
                      r2d = "Deviance Explained",
                      auc = "Area Under Curve",
                      mae = "Mean Absolute Error",
                      mce = "Mean Cross Entropy",
                      mse = "Mean Squared Error")
  )
  print(results_frame)
}

#' @export
print.predictive_gain <- function(x, digits = 3){
  results_frame <- as_data_frame(x[c("Model1", "Model2", "Delta")])
  results_frame$CI <- map_chr(
    x[[4]], ~ paste0("[", signif(.x[1],2), ", ", signif(.x[2],2), "]", sep = "")
  )
  results_frame <- as.data.frame(results_frame)
  metrics <- names(x$Model1)
  row.names(results_frame) <- map(
    metrics, ~ switch(.x,
                      rsq = "Variance Explained",
                      r2d = "Deviance Explained",
                      auc = "Area Under Curve",
                      mae = "Mean Absolute Error",
                      mce = "Mean Cross Entropy",
                      mse = "Mean Squared Error")
  )
  names(results_frame)[4] <- paste(names(x)[4], "for Delta")
  print(results_frame, digits = digits)
}

#' @export
print.summary_beset_elnet <- function(x, ...){
  cat("\n=======================================================",
      "\nBest Model:\n", sep = "")
  if(x$best$alpha < .25){
    cat("Primarily ridge ")
  } else if (x$best$alpha > .75){
    cat("Primarily lasso ")
  } else {
    cat("Mixture of ridge and lasso ")
  }
  cat("(alpha = ", x$best$alpha, " with lambda = ",
      x$best$best_lambda, ")", sep = "")
  coef_frame <- data_frame(
    variable = rownames(coef(x$best)),
    coef =  as.vector(coef(x$best, s = x$best$best_lambda))
  )
  if(!is.null(x$coef_ci)){
    coef_frame$conf.int <- as_data_frame(x$coef_ci) %>% transpose %>%
      simplify_all %>%
      map_chr(~paste0("[", signif(.x[1],3), ", ", signif(.x[2],3), "]",
                      collapse = ""))
  }
  coef_frame$stnd_coef <- NA
  coef_frame$stnd_coef[-1] <- coef_frame$coef[-1] * x$best$x_sd / x$best$y_sd
  coef_frame <- dplyr::filter(coef_frame, coef != 0)
  coef_frame <- dplyr::arrange(coef_frame, dplyr::desc(abs(stnd_coef)))
  coef_frame <- dplyr::mutate_if(coef_frame, is.numeric, ~ round(., 3))
  coef_frame <- as.data.frame(coef_frame)
  if(nrow(coef_frame) >= 1){
    cat("\n\nNon-zero coefficients ranked in order of importance:\n")
    print(coef_frame, quote = FALSE)
    cat("\n")
    cat(paste("Train-sample R-squared =", round(x$r2,2)))
    if(!is.null(x$r2_test)){
      cat(paste(", Test-sample R-squared =", round(x$r2_test, 2)))
    }
    cat("\n")
    cat(paste("Cross-validated R-squared = ", round(x$r2_cv$mean,2)))
  } else {
    cat("\n\nNo reliable predictors.")
  }
  cat("\n=======================================================")
}


#' @export
print.summary_beset_glm <- function(
  x, digits = max(3L, getOption("digits") - 3L, ...),
  signif.stars = getOption("show.signif.stars"), ...){
  object <- x
  cat("\n=======================================================",
      "\nBest Model:\n ", object$best_form, "\n")
  if(length(object$near_best) > 0){
    if(length(object$near_best) == 1){
      cat("\nNearly Equivalent Model:\n  ")
    } else {
      cat("\n", length(object$near_best), " Nearly Equivalent Models:\n  ",
          sep = "")
    }
    if(length(object$near_best) < 10){
      cat(object$near_best, sep = "\n  ")
    } else {
      cat(object$near_best[1:5], sep = "\n  ")
      cat("  ...\n   +", length(object$near_best) - 5, "more\n  ...\n")
    }
  }
  x <- object$best
  if (length(x$aliased) == 0L) {
    cat("\nNo Coefficients\n")
  } else {
    df <- if ("df" %in% names(x)) x[["df"]] else NULL
    if (!is.null(df) && (nsingular <- df[3L] - df[1L])){
      cat("\nCoefficients: (", nsingular,
          " not defined because of singularities)\n", sep = "")
    } else cat("\nCoefficients:\n")
    coefs <- signif(x$coefficients[, 1, drop = FALSE], digits)
    print(coefs)
  }
  cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ",
      format(x$dispersion), ")\n\n", sep = "")
  if (nzchar(mess <- stats::naprint(x$na.action)))
    cat("  (", mess, ")\n", sep = "")
  cat("Log-likelihood: ", formatC(x$loglik, digits = digits), " on ",
      x$loglik_df, " Df\nAIC: ",
      format(x$aic, digits = max(4L, digits + 1L)),
      "\n\n", "Number of Fisher Scoring iterations: ", x$iter,
      "\n", sep = "")
  correl <- object$best$correlation
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
  cat("\n")
  cat(paste("Train-sample R-squared =", round(object$r2,2)))
  if(!is.null(object$r2_test)){
    cat(paste(", Test-sample R-squared =", round(object$r2_test, 2)))
  }
  cat("\n")
  cat(paste("Cross-validated R-squared = ", round(object$r2_cv$mean,2)))
  cat("\n=======================================================")
}

#' @export
print.summary_nested_beset <- function(x, standardize = TRUE, metric = "rsq",
                                       ...){
  stnd <- if(standardize) "standardized" else "unstandardized"
  n_folds <- x$parameters$n_folds
  n_reps <- x$parameters$n_reps
  oneSE <- x$parameters$oneSE
  family <- x$parameters$family
  selection_metric <- x$parameters$metric
  selection_metric <- switch(x$parameters$metric,
                             auc = "Area Under Curve",
                             mae = "Mean Absolute Error",
                             mce = "Mean Cross Entropy",
                             mse = "Mean Squared Error")
  cat("\nResults of nested ", n_folds, "-fold cross-validation ", sep = "")
  if(n_reps > 1){
    cat("repeated ", n_reps, " times", sep = "")
  }
  cat("\n=======================================================\n")
  if(family == "negbin"){
    tune_frame <- data_frame(
      `Mean` = x$parameters$theta$mean,
      `S.E.` = x$parameters$theta$btwn_fold_se
      )
    if(n_reps > 1){
      theta_range <- x$parameters$theta$btwn_rep_range
      tune_frame$Range <- paste(signif(theta_range [1], 3),
                                signif(theta_range [2], 3), sep = " - ")
    }
    tune_frame <- tune_frame %>%
      mutate(Mean = signif(Mean, 3),
             S.E. = signif(S.E., 2))
    tune_frame <- as.data.frame(tune_frame)
    row.names(tune_frame) = "theta"
    print(tune_frame)
  }
  if(oneSE){
    cat("\nSimplest models within",
        "\n1 SE of best cross-validation ", selection_metric, ":\n", sep = "")
  } else {
    cat("\nModels with best cross-validation", selection_metric,
        ":\n")
  }
  form_frame <- as.data.frame(x$param$form)
  form_frame$Freq <- form_frame$Freq / sum(form_frame$Freq) * 100
  form_frame$Freq <- paste("(", round(form_frame$Freq), "%)", sep = "")
  names(form_frame)[1:2] <- ""
  print(form_frame, row.names = FALSE, width = 30)
  coef_frame <- data_frame(
    Coef. =  map_dbl(x$coefs[[stnd]], "mean"),
    S.E. = map_dbl(x$coefs[[stnd]], "btwn_fold_se"),
  )
  if(n_reps > 1){
    coef_frame$Min <- map_dbl(x$coefs[[stnd]], ~ .x$btwn_rep_range[1])
    coef_frame$Max <- map_dbl(x$coefs[[stnd]], ~ .x$btwn_rep_range[2])
  }
  coef_frame <- mutate_all(coef_frame, ~ round(., 3))
  coef_frame <- as.data.frame(coef_frame)
  row.names(coef_frame) = names(x$coefs[[stnd]])
  coef_frame <- coef_frame[coef_frame$Coef. != 0,]
  if(nrow(coef_frame) >= 1){
    if(standardize){
      names(coef_frame)[1] <- "Stnd.Coef."
      coef_frame <- coef_frame[order(abs(coef_frame$`Stnd.Coef`),
                                     decreasing = TRUE),]
    }
    cat("\n\nNon-zero coefficients")
    if(standardize) cat(" ranked in order of importance")
    cat(":\n")
    print(coef_frame, quote = FALSE)
  } else {
    cat("\n\nNo reliable predictors.")
  }
  cat("\n\nPrediction Metrics:\n")
  results_frame <- data_frame(
    `Mean` =  map_dbl(x$stats, ~ map_dbl(.x[metric], "mean")),
    `S.E.` = try(
      map_dbl(x$stats, ~ map_dbl(.x[metric], "btwn_fold_se")),
      silent = TRUE
    ))
  if(inherits(results_frame$S.E., "character")) results_frame$S.E. <- NA_real_
  if(n_reps > 1){
    results_frame$Min <- map_dbl(
      x$stats, ~ map_dbl(.x[metric], ~ .x$btwn_rep_range[1])
    )
    results_frame$Max <- map_dbl(
      x$stats, ~ map_dbl(.x[metric], ~ .x$btwn_rep_range[2])
    )
  }
  results_frame <- results_frame %>%
    mutate_all(~ signif(., 3)) %>%
    mutate_at(2, ~ signif(., 2))
  results_frame <- as.data.frame(results_frame)
  row.names(results_frame) <- c("Train Sample",
                                "CV-Tune Holdout",
                                "CV-Test Holdout")
  if(metric == "rsq" && family != "gaussian") metric <- "r2d"
  names(results_frame)[1] <- switch(metric,
                                    rsq = "Variance Explained",
                                    r2d = "Deviance Explained",
                                    auc = "Area Under Curve",
                                    mae = "Mean Absolute Error",
                                    mce = "Mean Cross Entropy",
                                    mse = "Mean Squared Error")
  print(results_frame)
  cat("=======================================================")
}

#' @export
print.summary_nested_elnet <- function(x, standardize = TRUE,
                                       metric = "rsq", ...){
  stnd <- if(standardize) "standardized" else "unstandardized"
  n_folds <- x$parameters$n_folds
  n_reps <- x$parameters$n_reps
  oneSE <- x$parameters$oneSE
  family <- x$parameters$family
  selection_metric <- x$parameters$metric
  selection_metric <- switch(x$parameters$metric,
                             auc = "Area Under Curve",
                             mae = "Mean Absolute Error",
                             mce = "Mean Cross Entropy",
                             mse = "Mean Squared Error")
  cat("\nResults of nested ", n_folds, "-fold cross-validation ", sep = "")
  if(n_reps > 1){
    cat("repeated ", n_reps, " times", sep = "")
  }
  cat("\n=======================================================")

  if(oneSE){
    cat("\nMost conservative tuning parameters within",
        "\n1 SE of best cross-validation ", selection_metric, ":\n", sep = "")
  } else {
    cat("\nTuning parameters with best cross-validation", selection_metric,
        ":\n")
  }
  tune_frame <- data_frame(
    `Mean` =  map_dbl(x$parameters[c("alpha", "lambda")], "mean"),
    `S.E.` = map_dbl(x$parameters[c("alpha", "lambda")], "btwn_fold_se")
  )
  tune_param <- c("alpha", "lambda")
  if(n_reps > 1){
    tune_frame$Range <- map_chr(
      x$parameters[tune_param],
      ~ paste(signif(.x$btwn_rep_range[1], 3),
              signif(.x$btwn_rep_range[2], 3), sep = " - "))
  }
  tune_frame <- dplyr::mutate_if(tune_frame, is.numeric, ~ signif(., 3))
  tune_frame <- as.data.frame(tune_frame)
  row.names(tune_frame) <- tune_param
  print(tune_frame)
  coef_frame <- data_frame(
    Coef. =  map_dbl(x$coefs[[stnd]], "mean"),
    S.E. = map_dbl(x$coefs[[stnd]], "btwn_fold_se")
  )
  if(n_reps > 1){
    coef_frame$Min <- map_dbl(x$coefs[[stnd]], ~ .x$btwn_rep_range[1])
    coef_frame$Max <- map_dbl(x$coefs[[stnd]], ~ .x$btwn_rep_range[2])
  }
  coef_frame <- mutate_all(coef_frame, ~ round(., 3))
  coef_frame <- as.data.frame(coef_frame)
  row.names(coef_frame) = names(x$coefs[[stnd]])
  coef_frame <- coef_frame[coef_frame$Coef. != 0,]
  if(nrow(coef_frame) >= 1){
    if(standardize){
      names(coef_frame)[1] <- "Stnd.Coef."
      coef_frame <- coef_frame[order(abs(coef_frame$`Stnd.Coef`),
                                     decreasing = TRUE),]
    }
    cat("\n\nNon-zero coefficients")
    if(standardize) cat(" ranked in order of importance")
    cat(":\n")
    print(coef_frame, quote = FALSE)
  } else {
    cat("\n\nNo reliable predictors.")
  }
  cat("\n\nPrediction Metrics:\n")
  results_frame <- data_frame(
    `Mean` =  map_dbl(x$stats, ~ map_dbl(.x[metric], "mean")),
    `S.E.` = map_dbl(x$stats, ~ map_dbl(.x[metric], "btwn_fold_se"))
  )
  if(n_reps > 1){
    results_frame$Min <- map_dbl(
      x$stats, ~ map_dbl(.x[metric], ~ .x$btwn_rep_range[1])
    )
    results_frame$Max <- map_dbl(
      x$stats, ~ map_dbl(.x[metric], ~ .x$btwn_rep_range[2])
    )
  }
  results_frame <- dplyr::mutate_all(results_frame, ~ signif(., 3))
  results_frame <- as.data.frame(results_frame)
  row.names(results_frame) <- c("Train Sample",
                                "CV-Tune Holdout",
                                "CV-Test Holdout")
  if(metric == "rsq" && family != "gaussian") metric <- "r2d"
  names(results_frame)[1] <- switch(metric,
                                    rsq = "Variance Explained",
                                    r2d = "Deviance Explained",
                                    auc = "Area Under Curve",
                                    mae = "Mean Absolute Error",
                                    mce = "Mean Cross Entropy",
                                    mse = "Mean Squared Error")
  print(results_frame)
  cat("=======================================================")
}

#' @export
print.variable_importance <- function(x, ...){
  print(plot(x))
}
