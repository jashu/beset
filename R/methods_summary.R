#' Summary Methods for \code{beset} Objects
#'
#' These functions are all methods for class \code{beset} objects.
#'
#' @param object An object for which a summary is desired.
#'
#' @param n_pred (Optional) \code{integer} number of predictors that the best
#' model should contain. If specified, all other arguments are ignored.
#'
#' @param metric \code{Character} string giving prediction metric on which to
#' base model selection. Can be one of \code{"auc"} for area under the (ROC)
#' curve (only available for binomial family), \code{"mae"} for mean absolute
#' error (not available for binomial family), \code{"mae"} for mean absolute
#' error, \code{"mce"} for mean cross entropy, or \code{"mse"} for mean squared
#' error. Default is \code{"auto"} which plots MSE for Gaussian-family models
#' and MCE for all other families.
#'
#' @param oneSE \code{Logical} indicating whether or not to use the "one
#' standard error" rule. If \code{TRUE} (default) the simplest model within one
#' standard error of the optimal model is returned. If \code{FALSE} the model
#' with the optimal cross-validation performance is returned.
#'
#' @param ... Additional arguments passed to model summary methods.
#'
#' @import purrr
#' @export
summary.beset <- function(object, n_pred = NULL, alpha = NULL, lambda = NULL,
                          metric = "auto", oneSE = TRUE, ...){
  if(inherits(object, "rf")){
    return(summary.beset_rf(object, ...))
  }
  metric <- tryCatch(
    match.arg(metric, c("auto", "aic", "auc", "mae", "mce", "mse", "rsq")),
    error = function(c){
      c$message <- gsub("arg", "metric", c$message)
      c$call <- NULL
      stop(c)
    }
  )
  tryCatch(
    if(
      (metric == "auc" && object$family != "binomial") ||
      (metric == "mae" && object$family == "binomial") ||
      (metric == "aic" && inherits(object, "elnet"))
    ) error = function(c){
      c$message <- paste(metric, "not available")
      c$call <- NULL
      stop(c)
    }
  )
  if(metric == "auto"){
    metric <- if(object$family == "gaussian") "mse" else "mce"
  }
  if(inherits(object, "elnet")){
    best_model <- get_best.beset_elnet(
      object, alpha = alpha, lambda = lambda, metric = metric[1],
      oneSE = oneSE)
    best_idx <- with(object$stats$fit,
                     which(alpha == best_model$alpha &
                             lambda == best_model$best_lambda))
    r2 <- object$stats$fit$rsq[best_idx]
    r2_test <- object$stats$test$rsq[best_idx]
    r2_cv <- purrr::flatten(object$stats$cv$rsq[best_idx])
    structure(list(best = best_model, r2 = r2, r2_cv = r2_cv, r2_test = r2_test),
              class = "summary_beset_elnet")
  } else {
    best_model <- beset:::get_best.beset_glm(
      object, n_pred = n_pred, metric = metric, oneSE = oneSE)
    loglik <- stats::logLik(best_model)
    best <- summary(best_model, ...)
    best$loglik <- loglik
    best$loglik_df <- attr(loglik, "df")
    best_form <- best_model$formula
    r2 <- object$stats$fit$rsq[object$stats$fit$form == best_form]
    r2_test <- if(!is.null(object$stats$test))
      object$stats$test$rsq[object$stats$test$form == best_form] else NULL
    r2_cv <- purrr::flatten(
      object$stats$cv$rsq[object$stats$cv$form == best_form])
    n_pred <- object$stats$fit$n_pred[object$stats$fit$form == best_form]
    near_equals <- object$stats$fit[object$stats$fit$n_pred == n_pred,]
    best_rsq <- max(near_equals$rsq, na.rm = T)
    near_best <- near_equals$form[near_equals$rsq > best_rsq - .01]
    near_best <- near_best[near_best != best_form]
    structure(list(best = best, best_form = best_form, near_best = near_best,
                   r2 = r2, r2_cv = r2_cv, r2_test = r2_test),
              class = "summary_beset_glm")
  }
}

#' @export
summary.nested <- function(object, metric = "auto", oneSE = TRUE, ...){
  metric <- tryCatch(
    match.arg(metric, c("auto", "aic", "auc", "mae", "mce", "mse", "rsq")),
    error = function(c){
      c$message <- gsub("arg", "metric", c$message)
      c$call <- NULL
      stop(c)
    }
  )
  tryCatch(
    if(
      (metric == "auc" && object$family != "binomial") ||
      (metric == "mae" && object$family == "binomial") ||
      (metric == "aic" && inherits(object, "elnet"))
    ) error = function(c){
      c$message <- paste(metric, "not available")
      c$call <- NULL
      stop(c)
    }
  )
  if(metric == "auto"){
    metric <- if(object$family == "gaussian") "mse" else "mce"
  }
  if(inherits(object, "elnet"))
    summary.nested_elnet(object, metric = metric, oneSE = oneSE, ...)
  else
    summary.nested_beset(object, metric = metric, oneSE = oneSE, ...)
}

summary.nested_elnet <- function(object, metric, oneSE, ...){
  family <- object$family
  n_folds <- object$n_folds
  n_reps <- object$n_reps
  repeats <- paste("Rep", 1:n_reps, "$", sep = "")
  rep_idx <- map(repeats, ~ grepl(.x, names(object$beset)))
  best_models <- map(
    object$beset,
    ~ beset:::get_best.beset_elnet(.x,  metric = metric, oneSE = oneSE)
  )
  alphas <- map_dbl(best_models, "alpha")
  lambdas <- map_dbl(best_models, "best_lambda")
  betas <- map2(best_models, lambdas, ~ as.vector(coef(.x, s = .y)))
  stnd <- map(best_models, ~.x$x_sd / .x$y_sd)
  stnd_betas <- map2(betas, stnd, ~ .x[-1] * .y) %>% transpose %>%
    simplify_all %>% as_tibble
  betas <- betas %>% transpose %>% simplify_all
  names(betas) <- best_models[[1]] %>% coef %>% row.names
  betas <- as_tibble(betas)
  betas <- list(
    mean = map(betas, mean),
    btwn_fold_se = map(betas, ~ sd(.x) / sqrt(n_folds)),
    btwn_rep_range = map(betas, function(x)
      map(rep_idx, ~ mean(x[.x])) %>% range)
  ) %>% transpose
  stnd_betas <- list(
    mean = map(stnd_betas, mean),
    btwn_fold_se = map(stnd_betas, ~ sd(.x) / sqrt(n_folds)),
    btwn_rep_range = map(stnd_betas, function(x)
      map(rep_idx, ~ mean(x[.x])) %>% range)
  ) %>% transpose
  best_idx <- pmap_int(
    list(s = map(object$beset, ~.x$stats$test), a = alphas, l = lambdas),
    function(s, a, l) with(s, which(alpha == a & lambda == l)))
  fit_stats <- map2_df(object$beset, best_idx, ~.x$stats$fit[.y,])
  fit_stats <- list(
    mean = map(fit_stats[3:6], ~mean(., na.rm = TRUE)),
    btwn_fold_se = map(fit_stats[3:6], ~ sd(.x, na.rm = TRUE) / sqrt(n_folds)),
    btwn_rep_range = map(fit_stats[, 3:6], function(x)
      map(rep_idx, ~ mean(x[.x])) %>% range)
  ) %>% transpose
  cv_stats <- map2_df(object$beset, best_idx, ~.x$stats$cv[.y,])
  cv_stats <- list(
    mean = map(cv_stats[3:6], ~ mean(map_dbl(., "mean"), na.rm = TRUE)),
    btwn_fold_se = map(cv_stats[3:6], ~ mean(map_dbl(., "btwn_fold_se"),
                                             na.rm = TRUE)),
    btwn_rep_range = map(
      cv_stats[, 3:6], function(x)
        map(rep_idx, ~ mean(map_dbl(x, "mean")[.x])) %>% range)
  ) %>% transpose
  validation_metrics <- validate(object, metric = metric, oneSE = oneSE)
  test_stats <- validation_metrics$stats
  best_lambda <- list(
    mean = map_dbl(rep_idx, ~ mean(lambdas[.x])) %>% mean,
    btwn_fold_se = sd(lambdas) / sqrt(n_folds),
    btwn_rep_range = map_dbl(rep_idx, ~ mean(lambdas[.x])) %>% range
  )
  best_alpha <- list(
    mean = map_dbl(rep_idx, ~ mean(alphas[.x])) %>% mean,
    btwn_fold_se = sd(alphas) / sqrt(n_folds),
    btwn_rep_range = map_dbl(rep_idx, ~ mean(alphas[.x])) %>% range
  )
  structure(list(
    stats = list(fit = fit_stats, cv = cv_stats, test = test_stats),
    parameters = c(list(alpha = best_alpha, lambda = best_lambda),
                   validation_metrics$parameters),
    coefs = list(unstandardized = betas, standardized = stnd_betas)
    ), class = "summary_nested_elnet"
    )
}

summary.nested_beset <- function(object, metric = "auto", oneSE = TRUE, ...){
  family <- object$family
  n_folds <- object$n_folds
  n_reps <- object$n_reps
  repeats <- paste("Rep", 1:n_reps, "$", sep = "")
  rep_idx <- purrr::map(repeats, ~ grepl(.x, names(object$fold_assignments)))
  best_models <- purrr::map(
    object$beset, ~ beset:::get_best.beset_glm(.x,  metric = metric,
                                               oneSE = oneSE)
  )
  if(family != "gaussian"){
    # Menard's standardization of y for logistic models
    var_logit_yhat <- map_dbl(best_models, ~ var(predict(.x)))
    rsq <- map_dbl(best_models, ~ r2d(.x)$R2fit)
    y_stnd <- var_logit_yhat / rsq
  } else {
    y_stnd <- map_dbl(object$beset, ~ sd(.x$parameters$y))
  }
  betas <- matrix(0, ncol = ncol(object$beset[[1]]$parameters$x),
                  nrow = length(object$beset))
  colnames(betas) <- colnames(object$beset[[1]]$parameters$x)
  for(i in 1:nrow(betas)){
    j <- coef(best_models[[i]])
    betas[i, names(j)] <- j
  }
  betas <- as_data_frame(betas)
  stnd <- purrr::map2(
    object$beset, y_stnd, ~.x$parameters$x %>% apply(2, sd) / .y
  ) %>% transpose %>% simplify_all
  stnd_betas <- map2(betas, stnd, ~ .x * .y) %>% as_data_frame
  betas <- list(
    mean = map(betas, mean),
    btwn_fold_se = map(betas, ~ sd(.x) / sqrt(n_folds)),
    btwn_rep_range = map(betas, function(x)
      map(rep_idx, ~ mean(x[.x])) %>% range)
  ) %>% transpose
  stnd_betas <- list(
    mean = map(stnd_betas, mean),
    btwn_fold_se = map(stnd_betas, ~ sd(.x) / sqrt(n_folds)),
    btwn_rep_range = map(stnd_betas, function(x)
      map(rep_idx, ~ mean(x[.x])) %>% range)
  ) %>% transpose
  best_form <- map_chr(best_models, "formula")
  best_idx <- match(best_form, object$beset[[1]]$stats$fit$form)
  fit_stats <- map2_df(object$beset, best_form,
                       ~ filter(.x$stats$fit, form == .y))
  stat_names <- intersect(names(fit_stats),
                          c("auc", "mae", "mce", "mse", "rsq"))
  fit_stats <- list(
    mean = map(fit_stats[stat_names], ~mean(., na.rm = TRUE)),
    btwn_fold_se = map(fit_stats[stat_names],
                       ~ sd(.x, na.rm = TRUE) / sqrt(n_folds)),
    btwn_rep_range = map(fit_stats[stat_names], function(x)
      map(rep_idx, ~ mean(x[.x])) %>% range)
  ) %>% transpose
  cv_stats <- map2_df(object$beset, best_form,
                      ~ filter(.x$stats$cv, form == .y))
  cv_stats <- list(
    mean = map(cv_stats[stat_names], ~ mean(map_dbl(., "mean"), na.rm = TRUE)),
    btwn_fold_se = map(cv_stats[stat_names], ~ mean(map_dbl(., "btwn_fold_se"),
                                                    na.rm = TRUE)),
    btwn_rep_range = map(
      cv_stats[, stat_names], function(x)
        map(rep_idx, ~ mean(map_dbl(x, "mean")[.x])) %>% range)
  ) %>% transpose
  validation_metrics <- validate(object, metric = metric, oneSE = oneSE)
  test_stats <- validation_metrics$stats
  form <- table(best_form)[order(table(best_form), decreasing = TRUE)]
  n_pred <- map2_int(object$beset, best_idx, ~.x$stats$fit$n_pred[.y])
  n_pred <- list(
    median = median(n_pred),
    btwn_fold_iqr = IQR(n_pred),
    btwn_rep_range = map(rep_idx, ~ median(n_pred[.x])) %>% range
  )
  structure(
    list(
      stats = list(fit = fit_stats, cv = cv_stats, test = test_stats),
      parameters = c(list(n_pred = n_pred, form = form),
                validation_metrics$parameters),
      coefs = list(unstandardized = betas, standardized = stnd_betas)
      ), class = "summary_nested_beset"
  )
}

#' @export
summary.beset_rf <- function(object, ...){
  type <- object$forests[[1]]$type
  ntree <- object$forests[[1]]$ntree
  mtry <- object$forests[[1]]$mtry
  n_folds <- object$stats$parameters$n_folds
  n_reps <- object$stats$parameters$n_reps
  repeats <- paste("Rep", 1:n_reps, "$", sep = "")
  rep_idx <- purrr::map(repeats, ~ grepl(.x, names(object$forests)))
  stats_oob <- if (type == "classification") {
    map_dbl(object$forests, ~.x$err.rate[ntree, "OOB"])
  } else {
    map_dbl(object$forests, ~ .x$rsq[ntree])
  }
  stats_test <- if (type == "classification") {
    map_dbl(object$forests, ~.x$test$err.rate[ntree, "Test"])
  } else {
    map_dbl(object$forests, ~ .x$test$rsq[ntree])
  }
  oob_stats <- list(
    mean = mean(stats_oob),
    btwn_fold_se = sd(stats_oob) / sqrt(n_folds),
    btwn_rep_range = map(rep_idx, ~ mean(stats_oob[.x])) %>% range)
  cv_stats <- list(
    mean = mean(stats_test),
    btwn_fold_se = sd(stats_test) / sqrt(n_folds),
    btwn_rep_range = map(rep_idx, ~ mean(stats_test[.x])) %>% range)
  validation_metrics <- validate(object, metric = metric, oneSE = oneSE)
  test_stats <- validation_metrics$stats
  varimp <- importance(object)
  structure(
    list(
      stats = list(oob = oob_stats, cv = cv_stats, test = test_stats),
      parameters = c(
        list(type = type, ntree = ntree, mtry = mtry),
        validation_metrics$parameters),
      vars = importance(object)
    ), class = "summary_beset_rf")
}
