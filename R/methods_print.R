#' @importFrom utils tail
#' @import purrr
#' @import dplyr

#' @export
print.beset <- function(x, ...) print(summary(x, ...))

#' @export
print.cross_valid <- function(x, digits = 2, ...){
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
  results_frame <- data.frame(
    Mean =  map_dbl(x$stats, "mean"),
    S.E. = map_dbl(x$stats, "btwn_fold_se")
  )
  if(x$parameters$n_reps > 1){
    results_frame$Min <- map_dbl(x$stats, ~ .x$btwn_rep_range[1])
    results_frame$Max <- map_dbl(x$stats, ~ .x$btwn_rep_range[2])
  }
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
  printCoefmat(as.matrix(results_frame), digits = digits,
               cs.ind = 1L:ncol(results_frame), tst.ind = NULL,
               P.values = FALSE, has.Pvalue = FALSE)
}

#' @export
print.prediction_metrics <- function(x, digits = 3, ...){
  results_frame <- data.frame(
    Metric =  map_dbl(x, ~ .) %>% signif(digits) %>% as.character
  )
  metrics <- names(x)
  names(results_frame) <- NULL
  if(attr(x, "family") != "gaussian") metrics[4] <- "r2d"
  row.names(results_frame) <- map(
    metrics, ~ switch(.x,
                      rsq = "Variance Explained",
                      r2d = "Deviance Explained",
                      auc = "Area Under Curve",
                      mae = "Mean Absolute Error",
                      mce = "Mean Cross Entropy",
                      mse = "Mean Squared Error")
  )
  print(results_frame, digits = digits)
}

#' @export
print.predictive_gain <- function(x, digits = 3, ...){
  results_frame <- as_tibble(x[c("Model1", "Model2", "Delta")])
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
  coef_frame <- tibble(
    variable = rownames(coef(x$best)),
    coef =  as.vector(coef(x$best, s = x$best$best_lambda))
  )
  if(!is.null(x$coef_ci)){
    coef_frame$conf.int <- as_tibble(x$coef_ci) %>% transpose %>%
      simplify_all %>%
      map_chr(~paste0("[", signif(.x[1],3), ", ", signif(.x[2],3), "]",
                      collapse = ""))
  }
  coef_frame$`stnd coef` <- NA
  coef_frame$`stnd coef`[-1] <- coef_frame$coef[-1] * x$best$x_sd / x$best$y_sd
  coef_frame <- dplyr::filter(coef_frame, coef != 0)
  coef_frame <- dplyr::arrange(coef_frame, dplyr::desc(abs(`stnd coef`)))
  coef_frame <- dplyr::mutate_if(coef_frame, is.numeric, ~ round(., 3))
  coef_frame <- as.data.frame(coef_frame)
  row.names(coef_frame) = coef_frame$variable
  coef_frame <- coef_frame[-1]
  if(nrow(coef_frame) >= 1){
    cat("\n\nNon-zero coefficients ranked in order of importance:\n")
    printCoefmat(coef_frame, digits = 3, quote = FALSE, has.Pvalue = FALSE,
                 cs.ind = 1:2, zap.ind = 1:2, tst.ind = NULL)
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
  out <- map2_chr(form_frame$best_form, form_frame$Freq, paste)
  cat(format(out, justify = "right"), sep = "\n")
  coef_frame <- data_frame(
    Coef. =  map_dbl(x$coefs[[stnd]], "mean"),
    S.E. = map_dbl(x$coefs[[stnd]], "btwn_fold_se"),
  )
  if(n_reps > 1){
    coef_frame$Min <- map_dbl(x$coefs[[stnd]], ~ .x$btwn_rep_range[1])
    coef_frame$Max <- map_dbl(x$coefs[[stnd]], ~ .x$btwn_rep_range[2])
  }
  # coef_frame <- mutate_all(coef_frame, ~ round(., 3))
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
    printCoefmat(coef_frame, digits = 3, quote = FALSE, has.Pvalue = FALSE,
                 cs.ind = 1:4, zap.ind = 1:4, tst.ind = NULL)
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
  printCoefmat(results_frame, digits = 3, quote = FALSE, has.Pvalue = FALSE,
               cs.ind = 1:4, zap.ind = 1:4, tst.ind = NULL)
  cat("=======================================================")
}

#' @export
print.summary_nested_elnet <- function(
  x, standardize = TRUE, metric = "rsq", ...
){
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
  tune_frame <- tibble(
    Mean =  map_dbl(x$parameters[c("alpha", "lambda")], "mean"),
    S.E. = map_dbl(x$parameters[c("alpha", "lambda")], "btwn_fold_se")
  )
  tune_param <- c("alpha", "lambda")
  if(n_reps > 1){
    tune_frame$Min <- map_dbl(x$parameters[tune_param], ~.x$btwn_rep_range[1])
    tune_frame$Max <- map_dbl(x$parameters[tune_param], ~.x$btwn_rep_range[2])
  }
  tune_frame <- dplyr::mutate_all(tune_frame, ~ round(., 3))
  tune_frame <- dplyr::mutate_all(tune_frame, ~ zapsmall(., 3))
  tune_frame <- as.data.frame(tune_frame)
  row.names(tune_frame) <- tune_param
  printCoefmat(tune_frame, P.values = FALSE, has.Pvalue = FALSE)
  coef_frame <- tibble(
    Coef. =  map_dbl(x$coefs[[stnd]], "mean"),
    S.E. = map_dbl(x$coefs[[stnd]], "btwn_fold_se")
  )
  if(n_reps > 1){
    coef_frame$Min <- map_dbl(x$coefs[[stnd]], ~ .x$btwn_rep_range[1])
    coef_frame$Max <- map_dbl(x$coefs[[stnd]], ~ .x$btwn_rep_range[2])
  }
  # coef_frame <- mutate_all(coef_frame, ~ round(., 3))
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
    printCoefmat(coef_frame, digits = 3, quote = FALSE, has.Pvalue = FALSE,
                 cs.ind = 1:4, zap.ind = 1:4, tst.ind = NULL)
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
  results_frame <- dplyr::mutate_all(results_frame, ~ round(., 3))
  results_frame <- dplyr::mutate_all(results_frame, ~ zapsmall(., 3))
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
  printCoefmat(results_frame, P.values = FALSE, has.Pvalue = FALSE)
  cat("=======================================================")
}

#' @export
print.variable_importance <- function(x, ...){
  print(plot(x))
}

#' @export
print.summary_beset_rf <- function(x, ...){
  n_folds <- x$parameters$n_folds
  n_reps <- x$parameters$n_reps
  family <- x$parameters$family
  type <- x$parameters$type
  ntree <- x$parameters$ntree
  mtry <- x$parameters$mtry


  cat("Type of random forest: ", type, "\n", sep = "")
  cat("Number of trees: ", ntree, "\n", sep = "")
  cat("No. of variables tried at each split: ", mtry, sep = "")
  cat("\n=======================================================\n")
  if (type == "classification") {
      cat("OOB estimate of error rate: ",
          round(x$stats$oob$mean * 100, digits = 2), "%\n", sep = "")
      cat("CV estimate of error rate: ",
          round(x$stats$cv$mean * 100, digits = 2), "%\n\n", sep = "")
  } else {
      cat("OOB estimate of % Var explained: ",
          round(100 * x$stats$oob$mean, digits = 2), "\n", sep = "")
      cat("CV estimate of % Var explained: ",
          round(100 * x$stats$cv$mean, digits = 2), "\n\n", sep = "")
  }
  var_imp <- x$vars %>% arrange(desc(importance))
  coef_frame <- data_frame(
    Importance =  var_imp$importance,
    Min = var_imp$min_import,
    Max = var_imp$max_import)
  coef_frame <- mutate_all(coef_frame, ~ round(., 3))
  coef_frame <- as.data.frame(coef_frame)
  row.names(coef_frame) <- var_imp$variable
  printCoefmat(coef_frame, digits = 3, quote = FALSE, has.Pvalue = FALSE,
               cs.ind = 1:3, zap.ind = 1:3, tst.ind = NULL)
  cat("\n\nPrediction Metrics\n")
  cat("(Results of ", n_folds, "-fold cross-validation ", sep = "")
  if(n_reps > 1){
    cat("repeated ", n_reps, " times", sep = "")
  }
  cat(")\n")
  results_frame <- results_frame <- data_frame(
    Mean =  map_dbl(x$stats$test, "mean"),
    S.E. = map_dbl(x$stats$test, "btwn_fold_se")
  )
  if(n_reps > 1){
    results_frame$Min <- map_dbl(x$stats$test, ~ .x$btwn_rep_range[1])
    results_frame$Max <- map_dbl(x$stats$test, ~ .x$btwn_rep_range[2])
  }
  results_frame <- dplyr::mutate_all(results_frame, ~ signif(., 3))
  results_frame <- as.data.frame(results_frame)
  metrics <- names(x$stats$test)
  if(family != "gaussian") metrics[4] <- "r2d"
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
  cat("=======================================================")
}

#' @export
print.beset_rf <- function(x, ...){
  print(summary(x, ...))
}

#' @export
print.R2 <- function(x, digits = 2){
  cat("Fit R-squared: ", formatC(x$R2fit, digits = digits))
  if(!is.null(x$R2new)){
    cat(",\tPrediction R-squared: ", formatC(x$R2new, digits = digits))
  }
  if(!is.null(x$R2cv)){
    cat(",\tCross-valid R-squared: ", formatC(x$R2cv, digits = digits))
  }
  cat("\n")
}
