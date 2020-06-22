#' Get best model from a "beset" class object
#'
#' \code{get_best} is a generic function used to obtain the best performing
#' predictive model from \code{\link{beset}} objects based on information
#' or cross-validation metrics.
#'
#' @param object A \code{\link{beset}} object.
#'
#' @param n_pred (Optional) \code{integer} number of predictors that the best
#' model should contain. If specified, all other arguments are ignored.
#'
#' @param metric \code{Character} string giving prediction metric on which to
#' base model selection. Can be one of \code{"mce"} (mean cross entropy--the
#' default), \code{"mse"} (mean squared error), \code{"aic"} (Akaike informatio
#' criterion, not applicable for elastic net), \code{"auc"} (area under the ROC
#' curve--only applicable if response is binomial), or \code{"mae"} (mean
#' absolute error--only applicable if response is numeric).
#'
#' @param oneSE \code{Logical} indicating whether or not to use the "one
#' standard error" rule. If \code{TRUE} (default) the simplest model within one
#' standard error of the optimal model is returned. If \code{FALSE} the model
#' with the optimal cross-validation performance is returned.
#'
#' @param alpha (Optional) \code{numeric} value to force selection of
#' elastic-net model with the given \code{alpha} parameter. If left \code{NULL},
#' the best value of \code{alpha} will be chosen using the cross-validation
#' \code{metric} and \code{oneSE} rule.
#'
#' @param lambda (Optional) \code{numeric} value to force selection of elastic-
#' net model with the given \code{lambda} parameter. If left \code{NULL},
#' the best value of \code{lambda} will be chosen using the cross-validation
#' \code{metric} and \code{oneSE} rule.

#' @export
get_best <- function(object, ...){
  UseMethod("get_best")
}


#' @export
get_best.glm <- function(object, n_pred = NULL, metric, oneSE = TRUE, ...){
  minimize <- TRUE; if(metric == "auc") minimize <- FALSE
  cv_stats <- object$stats$cv
  if(is.null(n_pred)){
    if(metric == "aic"){
      n_pred <- object$stats$fit$n_pred[which.min(object$stats$fit$aic)]
    } else {
      metric <- switch(
        metric,
        auc = cv_stats$auc,
        mae = cv_stats$mae,
        mce = cv_stats$mce,
        mse = cv_stats$mse
      )
      boundary <- find_boundary(metric, oneSE, minimize)
      n_pred <- get_best_par(par_name = "n_pred",
                             cv_stats = cv_stats,
                             metric = metric,
                             boundary = boundary,
                             minimize = minimize)
    }
  }
  best_idx <- which(cv_stats$n_pred == n_pred)
  fitter <- if(object$family == "negbin") "glm_nb" else "glm.fit"
  col_idx <- purrr::as_vector(cv_stats$pred_idx[best_idx])
  model_data <- c(
    list(x = object$parameters$x[, col_idx]),
    object$parameters[-1])
  best_model <- do.call(fitter, model_data)
  if(!inherits(best_model, "glm")) class(best_model) <- c("glm", "lm")
  best_model$formula <- cv_stats$form[best_idx]
  if(is.null(names(best_model$coefficients)))
    names(best_model$coefficients) <- "(Intercept)"
  best_model
}

#' @export
get_best.elnet <- function(
  object, metric, alpha = NULL, lambda = NULL, oneSE = TRUE, ...
){
  if(is.null(alpha) || is.null(lambda)){
    minimize <- TRUE; if(metric == "auc") minimize <- FALSE
    cv_stats <- object$stats$cv
    metric <- switch(metric,
                     auc = cv_stats$auc,
                     mae = cv_stats$mae,
                     mce = cv_stats$mce,
                     mse = cv_stats$mse)
    boundary <- find_boundary(metric, oneSE, minimize)
  }
  # if alpha value given, check to make sure this value was cross-validated
  best_alpha <- if(is.null(alpha)){
    get_best_par(
      par_name = "alpha", cv_stats = cv_stats, metric = metric,
      boundary = boundary, minimize = minimize
    )
  } else {
    if(!alpha %in% object$stats$fit$alpha)
      stop(paste(alpha, "was not among the cross-validated alpha values:\n",
                 paste0(unique(object$stats$fit$alpha), collapse = ", ")))
    alpha
  }
  best_lambda <- if(is.null(lambda)){
    idx <- which(cv_stats$alpha == best_alpha)
    get_best_par(
      par_name = "lambda", cv_stats = cv_stats[idx,], metric = metric[idx],
      boundary = boundary, minimize = minimize)
  } else {
  # if lambda value given, check to make sure this value was cross-validated
    valid_lambdas <- object$stats$fit %>% filter(alpha == best_alpha)
    lambda_range <- range(valid_lambdas$lambda)
    if(!between(lambda, lambda_range[1], lambda_range[2]))
      stop(
        paste(
          lambda, "outside the range of cross-validated lambda values:\n",
          paste0(lambda_range, collapse = " - ")
        )
      )
    valid_lambdas$lambda[which.min(abs(valid_lambdas$lambda - lambda))]
  }
  model_data <- object$parameters
  model_data$alpha <- best_alpha
  best_model <- do.call(glmnet::glmnet, model_data)
  best_model$alpha <- best_alpha
  best_model$best_lambda <- best_lambda
  x <- model_data$x
  y <- model_data$y
  if(is.factor(y)) y <- as.integer(y) - 1
  best_model$x_sd <- apply(x, 2, sd)
  best_model$y_sd <- if(object$family != "gaussian"){
    var_logit_yhat <- var(
      predict(best_model, newx = x, s = best_lambda, type = "link",
              newoffset = model_data$offset) %>% as.vector
    )
    rsq <- r2d(best_model)$R2fit
    var_logit_yhat / rsq
  } else sd(y)
  best_model$terms <- object$parameters$terms
  class(best_model) <- c("best", class(best_model))
  best_model
}

find_boundary <- function(error_metric, oneSE, minimize){
  best_error <- if(minimize){
    which.min(purrr::map_dbl(error_metric, "mean"))
  } else {
    which.max(purrr::map_dbl(error_metric, "mean"))
  }
  best_metric <- error_metric[[best_error]]
  out <- best_metric$mean
  if(oneSE){
    if(minimize){
      out <- out + best_metric$btwn_fold_se
    } else {
      out <- out - best_metric$btwn_fold_se
    }
  }
  out
}

get_best_par <- function(par_name, cv_stats, metric, boundary, minimize){
  best_1SE <- purrr::map_lgl(
    metric, ~ !is.na(.x$mean) &
      if(minimize){
        .x$mean <= boundary
      } else {
        .x$mean >= boundary
      }
    )
  if(par_name == "n_pred")
    min(cv_stats[[par_name]][best_1SE]) else max(cv_stats[[par_name]][best_1SE])
}
