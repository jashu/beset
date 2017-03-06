#' Summary Methods for \code{beset} Objects
#'
#' @inheritParams base::summary
#'
#' @param object An object of class \code{"beset_glm"},
#' \code{"beset_zeroinfl"}, or \code{"beset_elnet"}
#'
#' @param metric Prediction metric on which to based model selection. Can be one
#' of "MCE" (mean cross entropy), "MSE" (mean squared error), or "R2"
#' R-squared).
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
summary.beset_elnet <- function(object, metric = "MCE", oneSE = TRUE,
                                n_cores = 2, ...){
  metric <- match.arg(metric, c("MCE", "MSE", "R2"))
  if(is.na(metric)) stop("invalid 'metric' argument")
  best_metric <- switch(
    metric,
    MCE = min(object$stats$xval$MCE, na.rm = T)[1],
    MSE = min(object$stats$xval$MSE, na.rm = T)[1],
    R2 = max(object$stats$xval$R2, na.rm = T)[1]
  )
  best_idx <- switch(
    metric,
    MCE = which.min(object$stats$xval$MCE),
    MSE = which.min(object$stats$xval$MSE),
    R2 = which.max(object$stats$xval$R2)
  )
  best_alpha <- object$stats$xval$alpha[best_idx]
  best_lambda <- object$stats$xval$lambda[best_idx]
  if(oneSE){
    boundary <- switch(
      metric,
      MCE = best_metric + object$stats$xval$MCE_SE[best_idx],
      MSE = best_metric + object$stats$xval$MSE_SE[best_idx],
      R2 = best_metric - object$stats$xval$R2_SE[best_idx])
    best_1SE <- switch(
      metric,
      MCE = object$stats$xval[object$stats$xval$MCE < boundary,],
      MSE = object$stats$xval[object$stats$xval$MSE < boundary,],
      R2 = object$stats$xval[object$stats$xval$R2 > boundary,]
    )
    best_alpha <- max(best_1SE$alpha)
    best_1SE <- dplyr::filter(best_1SE, alpha == best_alpha)
    best_lambda <- max(best_1SE$lambda)
    best_idx <- which(object$stats$fit$alpha == best_alpha &
                        object$stats$fit$lambda == best_lambda)
  }
  best_model <- object$model_fits[[as.character(best_alpha)]]
  coefs <- coef(best_model, s = best_lambda)
  var_imp <- dplyr::data_frame(variable = rownames(coefs)[-1], coef = 0)
  var_imp$coef[coefs@i[-1]] <- coefs@x[-1]
  var_imp <- dplyr::arrange(var_imp, dplyr::desc(abs(coef)))
  R2 <- object$stats$fit$R2[best_idx]
  R2_test <- object$stats$test$R2[best_idx]
  x <- object$model_params$x
  y <- object$model_params$y
  family <- object$model_params$family
  cv_R2 <- sapply(object$xval_params$fold_ids, function(i){
    fit <- update(best_model, x = x[i,], y = y[i], family = family,
                  alpha = best_alpha)
    y_hat <- predict(fit, x[-i,], best_lambda, "response")
    predict_metrics_(y[-i], y_hat, family)$R_squared
  })
  R2_boot <- boot::boot(cv_R2, function(x,i) median(x[i]), 1000)
  boot_CI <- boot::boot.ci(R2_boot, type = "bca")
  R2_cv <- structure(
    list(cv_R2 = median(cv_R2, na.rm = TRUE), `95% CI` = boot_CI$bca[4:5],
         R2 = cv_R2), class = "cv_R2")

  structure(list(best = best_model, best_alpha = best_alpha,
                 best_lambda = best_lambda, var_imp = var_imp,
                 R2 = R2, R2_cv = R2_cv, R2_test = R2_test),
            class = "summary_beset_elnet")
}

#' @export
#' @rdname summary.beset
summary.beset_glm <- function(object, metric = "MCE", n_pred = NULL,
                              oneSE = TRUE, n_cores = 2, ...){
  metric <- tryCatch(match.arg(metric, c("AIC", "MCE", "MSE", "R2")),
                     error = function(c){
                       c$message <- gsub("arg", "metric", c$message)
                       c$call <- NULL
                       stop(c)
                     })
  best_model <- object$best_AIC
  best_form <- object$stats$fit$form[1]
  if(!is.null(n_pred)){
    best_form <- object$stats$cv$form[object$stats$cv$n_pred == n_pred]
    best_model <- if(class(best_model)[1] == "negbin"){
      update(best_model, best_form, data = object$model_data)
    } else {
      update(best_model, best_form, family = best_model$family,
             data = object$model_data)
    }
  }
  if(metric != "AIC" && is.null(n_pred)){
    best_metric <- switch(
      metric,
      MCE = min(object$stats$cv$MCE, na.rm = T)[1],
      MSE = min(object$stats$cv$MSE, na.rm = T)[1],
      R2 = max(object$stats$cv$R2, na.rm = T)[1]
    )
    best_idx <- switch(
      metric,
      MCE = which.min(object$stats$cv$MCE),
      MSE = which.min(object$stats$cv$MSE),
      R2 = which.max(object$stats$cv$R2)
    )
    best_form <- object$stats$cv$form[best_idx]
    if(oneSE){
      boundary <- switch(
        metric,
        MCE = best_metric + object$stats$cv$MCE_SE[best_idx],
        MSE = best_metric + object$stats$cv$MSE_SE[best_idx],
        R2 = best_metric - object$stats$cv$R2_SE[best_idx])
      best_1SE <- switch(
        metric,
        MCE = object$stats$cv[object$stats$cv$MCE < boundary,],
        MSE = object$stats$cv[object$stats$cv$MSE < boundary,],
        R2 = object$stats$cv[object$stats$cv$R2 > boundary,]
      )
      best_form <- best_1SE$form[which.min(best_1SE$n_pred)]
    }
    best_model <- if(class(best_model)[1] == "negbin"){
      update(best_model, best_form, data = object$model_data)
    } else {
      update(best_model, best_form, family = best_model$family,
             data = object$model_data)
    }
  }
  p <- object$stats$fit$n_pred[object$stats$fit$form == best_form]
  loglik <- logLik(best_model)
  best <- c(summary(best_model), loglik = loglik,
            loglik_df = attr(loglik, "df"))
  R2 <- object$stats$fit$R2[object$stats$fit$form == best_form]
  R2_test <- object$stats$test$R2[object$stats$test$form == best_form]
  R2_cv <- do.call("cv_r2",
                   args = c(list(object = best_model, n_cores = n_cores),
                            object$cv_params))
  near_equals <- dplyr::filter(object$stats$fit, n_pred == p)
  if(metric == "AIC") metric <- "MCE"
  near_equals <- dplyr::select(near_equals, form,
                               metric = dplyr::starts_with(metric))
  if(metric == "MCE") near_equals$metric <- exp(-near_equals$metric)
  best_metric <- if(metric == "MSE") min(near_equals$metric, na.rm = T)[1] else
    max(near_equals$metric, na.rm = T)[1]
  boundary <- c(best_metric - .01*best_metric, best_metric + .01*best_metric)
  near_best <- near_equals$form[dplyr::between(
    near_equals$metric, boundary[1], boundary[2])]
  near_best <- near_best[near_best != best_form]

  structure(list(best = best, best_form = best_form, near_best = near_best,
                 R2 = R2, R2_cv = R2_cv, R2_test = R2_test),
            class = "summary_beset_glm")
}

#' @export
#' @rdname summary.beset
summary.beset_zeroinfl <- function(object, metric = "MCE", n_count_pred = NULL,
                                   n_zero_pred = NULL, oneSE = TRUE,
                                   n_cores = 2, ...){
  metric <- tryCatch(match.arg(metric, c("AIC", "MCE", "MSE", "R2")),
                     error = function(c){
                       c$message <- gsub("arg", "metric", c$message)
                       c$call <- NULL
                       stop(c)
                     })
  best_model <- object$best_AIC
  best_form <- object$stats$fit$form[1]
  if(!is.null(n_count_pred) && !is.null(n_zero_pred)){
    best_idx <- object$stats$cv$n_count_pred == n_count_pred &
      object$stats$cv$n_zero_pred == n_zero_pred
    best_form <- object$stats$cv$form[best_idx]
    best_model <- pscl::zeroinfl(formula = formula(best_form),
                                 data = object$model_data,
                                 weights = best_model$weights,
                                 offset = best_model$offset$count,
                                 dist = best_model$dist,
                                 link = best_model$link,
                                 control = best_model$control)
  }
  if(metric != "AIC" && is.null(n_zero_pred) && is.null(n_count_pred)){
    best_metric <- switch(
      metric,
      MCE = min(object$stats$cv$MCE, na.rm = T)[1],
      MSE = min(object$stats$cv$MSE, na.rm = T)[1],
      R2 = max(object$stats$cv$R2, na.rm = T)[1]
    )
    best_idx <- switch(
      metric,
      MCE = which.min(object$stats$cv$MCE),
      MSE = which.min(object$stats$cv$MSE),
      R2 = which.max(object$stats$cv$R2)
    )
    best_form <- object$stats$cv$form[best_idx]
    if(oneSE){
      boundary <- switch(
        metric,
        MCE = best_metric + object$stats$cv$MCE_SE[best_idx],
        MSE = best_metric + object$stats$cv$MSE_SE[best_idx],
        R2 = best_metric - object$stats$cv$R2_SE[best_idx])
      best_1SE <- switch(
        metric,
        MCE = object$stats$cv[object$stats$cv$MCE < boundary,],
        MSE = object$stats$cv[object$stats$cv$MSE < boundary,],
        R2 = object$stats$cv[object$stats$cv$R2 > boundary,]
      )
      best_1SE <- dplyr::filter(best_1SE,
                                (n_zero_pred + n_count_pred) ==
                                  min(n_zero_pred + n_count_pred))
      best_idx <- switch(
        metric,
        MCE = which.min(best_1SE$MCE),
        MSE = which.min(best_1SE$MSE),
        R2 = which.max(best_1SE$R2)
      )
      best_form <- best_1SE$form[best_idx]
    }
    best_model <- pscl::zeroinfl(formula = formula(best_form),
                                 data = object$model_data,
                                 weights = best_model$weights,
                                 offset = best_model$offset$count,
                                 dist = best_model$dist,
                                 link = best_model$link,
                                 control = best_model$control)
  }
  pz <- object$stats$fit$n_zero_pred[object$stats$fit$form == best_form]
  pc <- object$stats$fit$n_count_pred[object$stats$fit$form == best_form]
  loglik <- logLik(best_model)
  best <- c(summary(best_model), loglik = loglik,
            loglik_df = attr(loglik, "df"), aic = AIC(best_model))
  R2 <- object$stats$fit$R2[object$stats$fit$form == best_form]
  R2_test <- object$stats$test$R2[object$stats$test$form == best_form]
  R2_cv <- do.call("cv_r2",
                   args = c(list(object = best_model, n_cores = n_cores),
                            object$cv_params))
  near_equals <- dplyr::filter(object$stats$fit,
                               (n_zero_pred + n_count_pred) == (pc + pz))
  if(metric == "AIC") metric <- "MCE"
  near_equals <- dplyr::select(near_equals, form,
                               metric = dplyr::starts_with(metric))
  if(metric == "MCE") near_equals$metric <- exp(-near_equals$metric)
  best_metric <- if(metric == "MSE") min(near_equals$metric, na.rm = T)[1] else
    max(near_equals$metric, na.rm = T)[1]
  boundary <- c(best_metric - .01*best_metric, best_metric + .01*best_metric)
  near_best <- near_equals$form[dplyr::between(
    near_equals$metric, boundary[1], boundary[2])]
  near_best <- near_best[near_best != best_form]

  structure(list(best = best, best_form = best_form, near_best = near_best,
                 R2 = R2, R2_cv = R2_cv, R2_test = R2_test),
            class = "summary_beset_zeroinfl")
}
