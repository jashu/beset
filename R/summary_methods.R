#' Summary Methods for \code{beset} Objects
#'
#' @inheritParams base::summary
#'
#' @param object An object of class \code{"beset_glm"},
#' \code{"beset_zeroinfl"}, or \code{"beset_elnet"}
#'
#' @param metric Prediction metric on which to based model selection. Can be one
#' of \code{"aic"} (Akaike information criterion), \code{"mae"} (mean absolute
#' error), \code{"mce"} (mean cross entropy), \code{"mse"} (mean squared error),
#' or \code{"r2"} (R-squared). If not specified, \code{"mse"} will be used for
#' Gaussian-family models and \code{"mce"} will be used for all other families.
#'
#' @param n_pred Optional: number of predictors that the best model should
#' contain
#'
#' @param n_count_pred Optional: number of predictors that the count component
#' of a zero-inflation model should contain.
#'
#' @param n_zero_pred Optional: number of predictors that the zero component of
#' a zero-inflation model should contain.
#'
#' @param oneSE Logical flag indicating whether or not the one-standard-error
#' rule should be applied when determining the best model from the model metric
#'
#' @param n_cores Number of CPUs to use when computing a bootstrapped confidence
#' interval for the cross-validated estimate of R-squared
#'
#' @name summary.beset
NULL

#' @export
#' @rdname summary.beset
summary.beset_elnet <- function(object, metric = "mce", oneSE = TRUE,
                                n_cores = 2, ...){
  metric <- match.arg(metric, c("mce", "mse", "r2"))
  if(is.na(metric)) stop("invalid 'metric' argument")
  best_metric <- switch(
    metric,
    mce = min(object$stats$cv$mce, na.rm = T)[1],
    mse = min(object$stats$cv$mse, na.rm = T)[1],
    r2 = max(object$stats$cv$r2, na.rm = T)[1]
  )
  best_idx <- switch(
    metric,
    mce = which.min(object$stats$cv$mce),
    mse = which.min(object$stats$cv$mse),
    r2 = which.max(object$stats$cv$r2)
  )
  best_alpha <- object$stats$cv$alpha[best_idx]
  best_lambda <- object$stats$cv$lambda[best_idx]
  if(oneSE){
    boundary <- switch(
      metric,
      mce = best_metric + object$stats$cv$mce_SE[best_idx],
      mse = best_metric + object$stats$cv$mse_SE[best_idx],
      r2 = best_metric - object$stats$cv$r2_SE[best_idx])
    best_1SE <- switch(
      metric,
      mce = object$stats$cv[object$stats$cv$mce < boundary,],
      mse = object$stats$cv[object$stats$cv$mse < boundary,],
      r2 = object$stats$cv[object$stats$cv$r2 > boundary,]
    )
    best_alpha <- max(best_1SE$alpha)
    best_1SE <- dplyr::filter(best_1SE, alpha == best_alpha)
    best_lambda <- max(best_1SE$lambda)
    best_idx <- which(object$stats$fit$alpha == best_alpha &
                        object$stats$fit$lambda == best_lambda)
  }
  best_model <- object$model_fits[[as.character(best_alpha)]]
  coefs <- glmnet::coef.glmnet(best_model, s = best_lambda)
  var_imp <- dplyr::data_frame(variable = rownames(coefs)[-1], stand.coef = 0)
  var_imp$stand.coef[coefs@i[-1]] <- coefs@x[-1]
  var_imp <- dplyr::arrange(var_imp, dplyr::desc(abs(stand.coef)))
  r2 <- object$stats$fit$r2[best_idx]
  r2_test <- object$stats$test$r2[best_idx]
  x <- object$model_params$x
  y <- object$model_params$y
  if(is.factor(y)) y <- as.integer(y) - 1
  family <- object$model_params$family
  cv_r2 <- sapply(object$cv_params$fold_ids, function(i){
    fit <- stats::update(best_model, x = x[i,], y = y[i], family = family,
                  alpha = best_alpha)
    y_hat <- stats::predict(fit, x[-i,], best_lambda, "response")
    predict_metrics_(y[-i], y_hat, family)$R_squared
  })
  r2_boot <- boot::boot(cv_r2, function(x,i) mean(x[i]), 1000)
  boot_CI <- boot::boot.ci(r2_boot, type = "bca")
  r2_cv <- structure(
    list(cv_R2 = mean(cv_r2, na.rm = TRUE), `95% CI` = boot_CI$bca[4:5],
         r2 = cv_r2), class = "cv_R2")

  structure(list(best = best_model, best_alpha = best_alpha,
                 best_lambda = best_lambda, var_imp = var_imp,
                 r2 = r2, r2_cv = r2_cv, r2_test = r2_test),
            class = "summary_beset_elnet")
}

#' @export
#' @rdname summary.beset
summary.beset_glm <- function(object, metric = NULL, n_pred = NULL,
                              oneSE = TRUE, n_cores = 2, ...){
  if(is.null(metric)){
    metric <- if(object$params$family == "gaussian") "mse" else "mce"
  }
  metric <- tryCatch(match.arg(metric, c("aic", "mae", "mce", "mse", "r2")),
                     error = function(c){
                       c$message <- gsub("arg", "metric", c$message)
                       c$call <- NULL
                       stop(c)
                     })
  best_form <- object$stats$fit$form[which.min(object$stats$fit$aic)]
  if(!is.null(n_pred)){
    best_form <- object$stats$cv$form[object$stats$cv$n_pred == n_pred]
  } else if (metric != "aic"){
    best_idx <- switch(
      metric,
      mae = which.min(object$stats$cv$mae),
      mce = which.min(object$stats$cv$mce),
      mse = which.min(object$stats$cv$mse),
      r2 = which.max(object$stats$cv$r2)
    )
    best_form <- object$stats$cv$form[best_idx]
    best_metric <- switch(
      metric,
      mae = object$stats$cv$mae[best_idx],
      mce = object$stats$cv$mce[best_idx],
      mse = object$stats$cv$mse[best_idx],
      r2 = object$stats$cv$r2[best_idx]
    )
    if(oneSE){
      boundary <- switch(
        metric,
        mae = best_metric + object$stats$cv$mae_SE[best_idx],
        mce = best_metric + object$stats$cv$mce_SE[best_idx],
        mse = best_metric + object$stats$cv$mse_SE[best_idx],
        r2 = best_metric - object$stats$cv$r2_SE[best_idx])
      best_1SE <- switch(
        metric,
        mae = object$stats$cv[object$stats$cv$mae < boundary,],
        mce = object$stats$cv[object$stats$cv$mce < boundary,],
        mse = object$stats$cv[object$stats$cv$mse < boundary,],
        r2 = object$stats$cv[object$stats$cv$r2 > boundary,]
      )
      best_form <- best_1SE$form[which.min(best_1SE$n_pred)]
    }
  }
  best_model <- fit_glm(mf = object$model_data, form = best_form,
                        family = object$params$family,
                        link = object$params$link)
  p <- object$stats$fit$n_pred[object$stats$fit$form == best_form]
  loglik <- stats::logLik(best_model)
  best <- c(summary(best_model), loglik = loglik,
            loglik_df = attr(loglik, "df"))
  r2 <- object$stats$fit$r2[object$stats$fit$form == best_form]
  r2_test <- object$stats$test$r2[object$stats$test$form == best_form]
  r2_cv <- cv_r2(best_model, n_cores = n_cores,
                 n_folds = object$params$n_folds,
                 n_repeats = object$params$n_repeats,
                 seed = object$params$seed)
  near_equals <- dplyr::filter(object$stats$fit, n_pred == p)
  near_equals <- dplyr::arrange(near_equals, dplyr::desc(r2))
  best_metric <- max(near_equals$r2, na.rm = T)[1]
  boundary <- c(best_metric - .01)
  near_best <- near_equals$form[near_equals$r2 > boundary]
  near_best <- near_best[near_best != best_form]
  structure(list(best = best, best_form = best_form, near_best = near_best,
                 r2 = r2, r2_cv = r2_cv, r2_test = r2_test),
            class = "summary_beset_glm")
}

#' @export
#' @rdname summary.beset
summary.beset_zeroinfl <- function(object, metric = "mce", n_count_pred = NULL,
                                   n_zero_pred = NULL, oneSE = TRUE,
                                   n_cores = 2, ...){
  metric <- tryCatch(match.arg(metric, c("aic", "mae", "mce", "mse", "r2")),
                     error = function(c){
                       c$message <- gsub("arg", "metric", c$message)
                       c$call <- NULL
                       stop(c)
                     })
  best_form <- object$stats$fit$form[which.min(object$stats$fit$aic)]
  if(!is.null(n_count_pred) && !is.null(n_zero_pred)){
    best_idx <- object$stats$cv$n_count_pred == n_count_pred &
      object$stats$cv$n_zero_pred == n_zero_pred
    best_form <- object$stats$cv$form[best_idx]
  } else if (metric != "aic"){
    best_metric <- switch(
      metric,
      mae = min(object$stats$cv$mae, na.rm = T)[1],
      mce = min(object$stats$cv$mce, na.rm = T)[1],
      mse = min(object$stats$cv$mse, na.rm = T)[1],
      r2 = max(object$stats$cv$r2, na.rm = T)[1]
    )
    best_idx <- switch(
      metric,
      mae = which.min(object$stats$cv$mae),
      mce = which.min(object$stats$cv$mce),
      mse = which.min(object$stats$cv$mse),
      r2 = which.max(object$stats$cv$r2)
    )
    best_form <- object$stats$cv$form[best_idx]
    if(oneSE){
      boundary <- switch(
        metric,
        mae = best_metric + object$stats$cv$mae_SE[best_idx],
        mce = best_metric + object$stats$cv$mce_SE[best_idx],
        mse = best_metric + object$stats$cv$mse_SE[best_idx],
        r2 = best_metric - object$stats$cv$r2_SE[best_idx])
      best_1SE <- switch(
        metric,
        mae = object$stats$cv[object$stats$cv$mae < boundary,],
        mce = object$stats$cv[object$stats$cv$mce < boundary,],
        mse = object$stats$cv[object$stats$cv$mse < boundary,],
        r2 = object$stats$cv[object$stats$cv$r2 > boundary,]
      )
      best_1SE <- dplyr::filter(best_1SE,
                                (n_zero_pred + n_count_pred) ==
                                  min(n_zero_pred + n_count_pred))
      best_idx <- switch(
        metric,
        mae = which.min(best_1SE$mae),
        mce = which.min(best_1SE$mce),
        mse = which.min(best_1SE$mse),
        r2 = which.max(best_1SE$r2)
      )
      best_form <- best_1SE$form[best_idx]
    }
  }
  best_model <- pscl::zeroinfl(formula = stats::formula(best_form),
                               data = object$model_data,
                               dist = object$params$family,
                               link = object$params$link)
  pz <- object$stats$fit$n_zero_pred[object$stats$fit$form == best_form]
  pc <- object$stats$fit$n_count_pred[object$stats$fit$form == best_form]
  loglik <- stats::logLik(best_model)
  best <- c(summary(best_model), loglik = loglik,
            loglik_df = attr(loglik, "df"), aic = stats::AIC(best_model))
  r2 <- object$stats$fit$r2[object$stats$fit$form == best_form]
  r2_test <- object$stats$test$r2[object$stats$test$form == best_form]
  r2_cv <- cv_r2(best_model, n_cores = n_cores,
                 n_folds = object$params$n_folds,
                 n_repeats = object$params$n_repeats,
                 seed = object$params$seed)
  near_equals <- if(is.null(n_zero_pred) || is.null(n_count_pred)){
    dplyr::filter(object$stats$fit, (n_zero_pred + n_count_pred) == (pc + pz))
  } else {
    dplyr::filter(object$stats$fit, n_zero_pred == pz & n_count_pred == pc)
  }
  best_metric <- max(near_equals$r2, na.rm = T)[1]
  boundary <- best_metric - .01
  near_best <- near_equals$form[near_equals$r2 > boundary]
  near_best <- near_best[near_best != best_form]

  structure(list(best = best, best_form = best_form, near_best = near_best,
                 r2 = r2, r2_cv = r2_cv, r2_test = r2_test),
            class = "summary_beset_zeroinfl")
}
