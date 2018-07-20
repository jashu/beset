#' Plot Methods for \code{beset} Objects
#'
#' @param x An object of class \code{"beset_glm"},
#' \code{"beset_zeroinfl"}, or \code{"beset_elnet"}
#'
#' @param metric \code{Character} string giving the performance metric to plot.
#' Can be one of \code{"auc"} for area under the (ROC) curve (only available for
#' binomial family), \code{"mae"} for mean absolute error (not available for
#' binomial family), \code{"mce"} for mean cross entropy, \code{"mse"} for
#' mean squared error, or \code{"rsq"} for R-squared. Default is \code{"auto"}
#' which plots MSE for Gaussian-family models and MCE for all other families.
#'
#' @param se \code{Logical} indicating whether or not to plot error bars.
#' \code{TRUE} by default.
#' @import purrr
#' @import dplyr
#' @import ggplot2
#' @export
plot.beset <- function(x, metric = "auto", se = TRUE, ...){
  metric <- tryCatch(
    match.arg(metric, c("auto", "auc", "mae", "mce", "mse", "rsq")),
    error = function(c){
      c$message <- gsub("arg", "metric", c$message)
      c$call <- NULL
      stop(c)
    }
  )
  tryCatch(
    if(
      (metric == "auc" && x$family != "binomial") ||
      (metric == "mae" && x$family == "binomial")
    ) error = function(c){
      c$message <- paste(metric, "not available for", x$family, "models")
      c$call <- NULL
      stop(c)
    }
  )
  if(metric == "auto"){
    metric <- if(x$family == "gaussian") "mse" else "mce"
  }
  if(inherits(x, "nested")){
    data <- beset:::aggregate.nested(x, metric)
  } else {
    parameters <- intersect(
      c("lambda", "alpha", "n_pred", "form", metric), names(x$stats$fit)
    )
    cv <- x$stats$cv[parameters]
    cv$cv <- map_dbl(cv[[metric]], "mean")
    cv$cv_lower <- cv$cv - purrr::map_dbl(cv[[metric]], "btwn_fold_se")
    cv$cv_upper <- cv$cv + purrr::map_dbl(cv[[metric]], "btwn_fold_se")
    cv[[metric]] <- NULL
    train <- x$stats$fit[parameters]
    names(train)[3] <- "train"
    data <- inner_join(train, cv, by = names(train)[1:2])
    test <- x$stats$test[parameters]
    if(!is.null(test)){
      names(test)[3] <- "test"
      data <- inner_join(data, test, by = names(train)[1:2])
    }
  }
  if(any(is.na(c(data$cv_lower, data$cv_upper)))){
    se <- FALSE
    warning(
      paste("Error bars could not be plotted: ",
            "Standard errors not available for ", metric, ".", sep = "")
    )
  }
  names(data)[1] <- "x"
  if("alpha" %in% names(data)){
    data <- data %>% group_by(alpha) %>%
      filter(between(x, quantile(x, .1), quantile(x, .9))) %>% ungroup()
  }
  best_idx <- ifelse(grepl("^m", metric),
                     which.min(data$cv), which.max(data$cv))
  x_int <- data[[1]][best_idx]
  y_int <- ifelse(grepl("^m", metric),
                  data$cv_upper[best_idx], data$cv_lower[best_idx])
  if(metric == "rsq" && x$family != "gaussian") metric <- "r2d"
  x_lab <- if(inherits(x, "elnet")){
    xlab("Regularization parameter")
  } else {
    xlab("Number of Predictors")
  }
  y_lab <- switch(metric,
                  auc = ylab("Area Under the ROC curve"),
                  mae = ylab("Mean Absolute Error"),
                  mce = ylab("Mean Cross Entropy"),
                  mse = ylab("Mean Squared Error"),
                  rsq = ylab(bquote(~R^2)),
                  r2d = ylab(bquote(~R[D]^2)))
  color_legend <- c(
    "Train Sample" = "grey", "CV Holdout" = "red", "Test Holdout" = "blue"
  )
  scale_x <- if(inherits(x, "elnet")){
    scale_x_log10()
  } else {
    scale_x_continuous(breaks = 0:max(data$x))
  }
  p <- ggplot(data = data) + theme_bw() + x_lab + y_lab + scale_x +
    scale_color_manual(name = "", values = color_legend) +
    geom_line(aes(x = x, y = train, color = "Train Sample"), size = 1) +
    geom_line(aes(x = x, y = cv, color = "CV Holdout"), size = 1)
  if("test" %in% names(data)){
    p <- p + geom_line(aes(x = x, y = test, color = "Test Holdout"), size = 1)
  }
  if(se){
    p <- p +
      geom_hline(yintercept = y_int, linetype = "dashed") +
      geom_vline(xintercept = x_int, linetype = "dashed")
    if(inherits(x, "elnet")){
      p <- p + geom_ribbon(aes(x = x, ymin = cv_lower, ymax = cv_upper),
                           alpha = .25, color = NA, fill = "grey70")
    } else if(inherits(x, "nested")) {
     p <- p +
       geom_ribbon(aes(x = x, ymin = cv_lower, ymax = cv_upper),
                   alpha = .25, color = NA, fill = "grey70") +
       geom_ribbon(aes(x = x, ymin = train_lower, ymax = train_upper),
                   alpha = .25, color = NA, fill = "grey70") +
       geom_ribbon(aes(x = x, ymin = test_lower, ymax = test_upper),
                   alpha = .25, color = NA, fill = "grey70")
    } else {
      p <- p + geom_errorbar(aes(x = x, ymin = cv_lower, ymax = cv_upper),
                             width = 0.2, color = "red")
    }
  }
  if("alpha" %in% names(data)){
    p <- p + facet_grid(alpha ~ .)
  }
  p
}
