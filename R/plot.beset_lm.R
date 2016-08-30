#' Plot Method for the \code{beset_lm} Class
#'
#' \code{plot.beset_lm} takes the output of a \code{\link{beset_lm}} object and creates
#' a line plot using \code{\link[ggplot2]{ggplot}}.
#'
#' \code{plot.beset_lm} produces a plot showing \eqn{R^2} as a function of the
#' number of model parameters from a cross-validated best-subset-selection
#' procedure, as implemented by \code{\link{beset_lm}}. Line graphs are plotted
#' for the training, cross-validated, and, if included, independent-test
#' \eqn{R^2}, which are labeled in the legend as "Train", "CV", and "Test",
#' respectively. If more than one choice of \code{k} was used for cross-
#' validation, separate panels are graphed for each.
#'
#' @param object An object of class \code{beset_lm}.
#'
#' @param SE Logical indicating whether or not standard-error bars should be
#' plotted.
#'
#' @param title Optional character vector giving plot title.
#'
#' @import ggplot2
#'
#' @export
plot.beset_lm <- function(object, SE = TRUE, title = ""){
  data <- object$R2
  data$R2_train_lower <- with(data, R2_train - R2_train_SE)
  data$R2_train_upper <- with(data, R2_train + R2_train_SE)
  data$R2_cv_lower <- with(data, R2_cv - R2_cv_SE)
  data$R2_cv_upper <- with(data, R2_cv + R2_cv_SE)
  data$R2_test_lower <- with(data, R2_test - R2_test_SE)
  data$R2_test_upper <- with(data, R2_test + R2_test_SE)
  ymin <- min(0, data$R2_cv_lower, data$R2_test_lower, na.rm = T)
  ymax <- max(data$R2_cv_upper, data$R2_test_upper, data$R2_train_upper, 1,
              na.rm = T)
  xmax <- max(data$n_preds)
  color_legend <- c("Train" = "grey", "CV" = "red", "Test" = "blue")
  p <- ggplot(data = data) +
    theme_bw() +
    ggtitle(title) +
    xlab("Number of Predictors") +
    ylab(expression(R^{2})) +
    scale_x_continuous(breaks = 1:xmax) +
    scale_y_continuous(limits = c(ymin, ymax)) +
    scale_color_manual(name = "", values = color_legend) +
    geom_line(aes(x = n_preds, y = R2_train, color = "Train")) +
    geom_line(aes(x = n_preds, y = R2_cv, color = "CV"))
  if(SE){
    p <- p +
      geom_errorbar(aes(x = n_preds, ymin = R2_train_lower,
                        ymax = R2_train_upper), width = 0.2, color = "grey") +
      geom_errorbar(aes(x = n_preds, ymin = R2_cv_lower, ymax = R2_cv_upper),
                    width = 0.2, color = "red")
  }
  if(!any(is.na(data$R2_test))){
    p <- p + geom_line(aes(x = n_preds, y = R2_test, color = "Test"))
    if(SE){
      p <- p + geom_errorbar(aes(x = n_preds, ymin = R2_test_lower,
                                 ymax = R2_test_upper), width = 0.2,
                             color = "blue")
    }
  }
  p
}

#' @export
plot.beset_glm <- function(object, SE = TRUE, title = ""){
  class(object) <- "beset_lm"
  plot(object, SE, title)
}
