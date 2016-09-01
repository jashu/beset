#' Plot Method for the \code{beset_zeroinfl} Class
#'
#' \code{plot.beset_zeroinfl} takes the output of a \code{\link{beset_zeroinfl}} object and creates
#' a line plot using \code{\link[ggplot2]{ggplot}}.
#'
#' \code{plot.beset_zeroinfl} produces a plot showing \eqn{R^2} as a function o
#' of the number of count-model parameters (x axis) and zero-model parameters
#' (individual panels) from a cross-validated  best-subset-selection procedure,
#' as implemented by \code{\link{beset_zeroinfl}}. Line graphs are plotted for
#' the training, cross-validated, and, if included, independent-test
#' \eqn{R^2}, which are labeled in the legend as "Train", "CV", and "Test",
#' respectively.
#'
#' @param object An object of class \code{beset_zeroinfl}.
#'
#' @param SE Logical indicating whether or not standard-error bars should be
#' plotted.
#'
#' @param title Optional character vector giving plot title.
#'
#' @import ggplot2
#'
#' @export
plot.beset_zeroinfl <- function(object, SE = TRUE, title = ""){
  data <- object$best_subsets
  if(any(is.na(data$test_CE))) data <- dplyr::select(data, -test_CE)
  data <- tidyr::gather(data, type, CE, train_CE:cv_CE)
  data$type <- factor(gsub("_CE", "", data$type))
  data$n_zero_pred <- factor(data$n_zero_pred)
  data$cv_CE_lower <- NA
  data$cv_CE_upper <- NA
  data$cv_CE_lower[data$type == "cv"] <-
    with(data[data$type == "cv",], CE - cv_CE_SE)
  data$cv_CE_upper[data$type == "cv"] <-
    with(data[data$type == "cv",], CE + cv_CE_SE)
  xmax <- max(data$n_count_pred)
  p <- ggplot(data = data,
              aes(x = n_count_pred,
                  y = CE,
                  color = type)) +
    facet_wrap(~ n_zero_pred) +
    ggtitle(title) +
    xlab("Number of Predictors in Count Model") +
    ylab("Cross-Entropy Error") +
    scale_x_continuous(breaks = 0:xmax) +
    geom_line() +
    theme_bw() +
    theme(legend.title = element_blank())
  if(SE){
    p <- p + geom_errorbar(data = data[data$type == "cv",],
                           aes(x = n_count_pred, ymin = cv_CE_lower,
                               ymax = cv_CE_upper), width = 0.2)
  }
  footnote <- paste("Numbers at top of each panel indicate",
                     "number of predictors in zero-inflation model.")
  grid::grid.newpage()
  g <- gridExtra::arrangeGrob(p, bottom = grid::textGrob(
    footnote, x = 0.5,  just = "centre",
    gp = grid::gpar(fontface = "italic", fontsize = 10)))
  grid::grid.draw(g)
}

#' @export
