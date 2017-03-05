#' Plot Methods for \code{beset} Objects
#'
#' @param object An object of class \code{"beset_glm"},
#' \code{"beset_zeroinfl"}, or \code{"beset_elnet"}
#'
#' @param type Type of error to plot. Can be one of \code{"train"}, \code{"cv"},
#' or \code{"test"}.
#'
#' @param metric Which error metric to plot. Can be one of \code{"MCE"} for
#' mean cross entropy, \code{"MSE"} for mean squared error, or \code{"R2"} for
#' R-squared.
#'
#' @name plot.beset
#'
#' @import ggplot2
NULL

#' @export
#' @rdname plot.beset
plot.beset_elnet <- function(object, type = "cv", metric = "MCE"){
  data <- switch(
    type,
    train = dplyr::select(object$stats$fit, alpha, lambda,
                          dplyr::starts_with(metric)),
    cv = dplyr::select(object$stats$cv, alpha, lambda,
                         dplyr::starts_with(metric)),
    test = tryCatch(dplyr::select(object$stats$test, alpha, lambda,
                                  dplyr::starts_with(metric)),
                    error = function(c){
                      c$message <- "test data not found"
                      stop(c)
                    })
  )
  data$alpha <- factor(data$alpha)
  if(type == "cv"){
    names(data)[3:4] <- c("error", "se")
    data$cv_lower <- with(data, error - se)
    data$cv_upper <- with(data, error + se)
  } else{
    names(data)[3] <- "error"
  }
  title <- ggtitle(switch(type,
                          train = "Fit to training data",
                          cv = "Cross-validation error",
                          test = "Prediction of new data"))
  y_lab <- switch(metric,
                  MCE = ylab("Mean Cross Entropy"),
                  MSE = ylab("Mean Squared Error"),
                  R2 = ylab(bquote(~R^2)))

  p <- ggplot(data = data, aes(x = lambda, y = error, color = alpha)) +
    theme_bw() + title + xlab("Regularization parameter") + y_lab +
    scale_x_log10() + geom_line() + labs(color = "Mixing parameter")
  if(type == "cv"){
    p <- p +
      geom_errorbar(aes(x = lambda, ymin = cv_lower, ymax = cv_upper),
                    width = 0.1)
  }
  p
}

#' @export
#' @rdname plot.beset
plot.beset_glm <- function(object, metric = "MCE"){
  train <- dplyr::select(object$stats$fit, n_pred, form,
                         dplyr::starts_with(metric))
  names(train)[3] <- "train"
  cv <- dplyr::select(object$stats$cv, n_pred, form,
                        dplyr::starts_with(metric))
  names(cv)[3:4] <- c("cv", "cv_se")
  test <- try(dplyr::select(object$stats$test, n_pred, form,
                            dplyr::starts_with(metric)), silent = TRUE)
  if(class(test)[1] == "try-error"){
    test <- NULL
  } else {
    names(test)[3] <- "test"
  }
  train <- dplyr::filter(train, form %in% cv$form)
  data <- dplyr::left_join(train, cv, by = c("n_pred", "form"))
  if(!is.null(test)) data <- dplyr::left_join(data, test,
                                              by = c("n_pred", "form"))
  data$cv_lower <- with(data, cv - cv_se)
  data$cv_upper <- with(data, cv + cv_se)
  xmax <- max(data$n_pred)
  color_legend <- c("Train" = "grey", "CV" = "red", "Test" = "blue")
  if(metric == "R2" && object$best_AIC$family$family != "gaussian")
    metric <- "R2D"
  y_lab <- switch(metric,
                  MCE = ylab("Mean Cross Entropy"),
                  MSE = ylab("Mean Squared Error"),
                  R2 = ylab(bquote(~R^2)),
                  R2D = ylab(bquote(~R[D]^2)))
  p <- ggplot(data = data) +
    theme_bw() +
    xlab("Number of Predictors") +
    y_lab +
    scale_x_continuous(breaks = 0:xmax) +
    scale_color_manual(name = "", values = color_legend) +
    geom_line(aes(x = n_pred, y = train, color = "Train")) +
    geom_line(aes(x = n_pred, y = cv, color = "CV")) +
    geom_errorbar(aes(x = n_pred, ymin = cv_lower, ymax = cv_upper),
                  width = 0.2, color = "red")
  if(!is.null(test)){
    p <- p + geom_line(aes(x = n_pred, y = test, color = "Test"))
  }
  p
}

#' @export
#' @rdname plot.beset
plot.beset_zeroinfl <- function(object, type = "cv", metric = "MCE"){
  data <- switch(
    type,
    train = dplyr::select(object$stats$fit, form, n_count_pred, n_zero_pred,
                          dplyr::starts_with(metric)),
    cv = dplyr::select(object$stats$cv, form, n_count_pred, n_zero_pred,
                         dplyr::starts_with(metric)),
    test = tryCatch(dplyr::select(object$stats$test, form, n_count_pred,
                                  n_zero_pred, dplyr::starts_with(metric)),
                    error = function(c){
                      c$message <- "test data not found"
                      stop(c)
                    })
  )
  data$n_zero_pred <- factor(data$n_zero_pred)
  if(type == "cv"){
    names(data)[4:5] <- c("error", "se")
    data$cv_lower <- with(data, error - se)
    data$cv_upper <- with(data, error + se)
  } else{
    names(data)[4] <- "error"
  }
  title <- ggtitle(switch(type,
                          train = "Fit to training data",
                          cv = "Cross-validation error",
                          test = "Prediction of new data"))
  y_lab <- switch(metric,
                  MCE = ylab("Mean Cross Entropy (-loglik / N)"),
                  MSE = ylab("Mean Squared Error"),
                  R2 = ylab(bquote(~R^2)))
  p <- ggplot(data = data,
              aes(x = n_count_pred, y = error, color = n_zero_pred)) +
    theme_bw() + title + xlab("Number of Count-Component Predictors") +
    y_lab + scale_x_continuous(breaks = 0:max(data$n_count_pred)) +
    geom_line() + labs(color = "Number of\nZero-Component\nPredictors")
  if(type == "cv"){
    p <- p +
      geom_errorbar(aes(x = n_count_pred, ymin = cv_lower, ymax = cv_upper),
                    width = 0.2)
  }
  p
}
