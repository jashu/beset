#' Plot Method for the \code{beset_glm} Class
#'
#' \code{plot.beset_glm} takes the output of a \code{\link{beset_glm}} object and creates
#' a line plot using \code{\link[ggplot2]{ggplot}}.
#'
#' \code{plot.beset_glm} produces a plot showing \eqn{R^2} as a function of the
#' number of model parameters from a cross-validated best-subset-selection
#' procedure, as implemented by \code{\link{beset_glm}}. Line graphs are plotted
#' for the training, cross-validated, and, if included, independent-test
#' \eqn{R^2}, which are labeled in the legend as "Train", "CV", and "Test",
#' respectively.
#'
#' @param object An object of class \code{beset_glm}.
#'
#' @param SE Logical indicating whether or not standard-error bars should be
#' plotted.
#'
#' @param title Optional character vector giving plot title.
#'
#' @import ggplot2
#'
#' @export
plot.beset_glm <- function(object, SE = TRUE, title = ""){
  data <- object$best_subsets
  data$R2_cv_lower <- with(data, cv_R2 - cv_R2_SE)
  data$R2_cv_upper <- with(data, cv_R2 + cv_R2_SE)
  ymin <- min(data$train_R2, data$R2_cv_lower, data$test_R2, na.rm = T)
  ymax <- max(data$train_R2, data$R2_cv_upper, data$test_R2, 1, na.rm = T)
  xmax <- max(data$n_pred)
  color_legend <- c("Train" = "grey", "CV" = "red", "Test" = "blue")
  p <- ggplot(data = data) +
    theme_bw() +
    ggtitle(title) +
    xlab("Number of Predictors") +
    ylab(expression(R^{2})) +
    scale_x_continuous(breaks = 0:xmax) +
    scale_y_continuous(limits = c(ymin, ymax)) +
    scale_color_manual(name = "", values = color_legend) +
    geom_line(aes(x = n_pred, y = train_R2, color = "Train")) +
    geom_line(aes(x = n_pred, y = cv_R2, color = "CV"))
  if(SE){
    p <- p +
      geom_errorbar(aes(x = n_pred, ymin = R2_cv_lower, ymax = R2_cv_upper),
                    width = 0.2, color = "red")
  }
  if(!any(is.na(data$test_R2))){
    p <- p + geom_line(aes(x = n_pred, y = test_R2, color = "Test"))
  }
  p
}
