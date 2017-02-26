#' Plot Method for the \code{beset_zeroinfl} Class
#'
#' \code{plot.beset_zeroinfl} takes the output of a \code{\link{beset_zeroinfl}} object and creates
#' a line plot using \code{\link[ggplot2]{ggplot}}.
#'
#' \code{plot.beset_zeroinfl} produces a plot showing cross-validated prediction
#' error as a function of the number of model parameters from a
#' \code{\link{beset_zeroinfl}} object. Line graphs are plotted for the
#' training, cross-validated, and, if included, independent-test error, which
#' are labeled in the legend as "Train", "CV", and "Test", respectively.
#'
#' @param object An object of class \code{beset_zeroinfl}.
#'
#' @param SE Logical indicating whether or not standard-error bars should be
#' plotted.
#'
#' @import ggplot2
#'
#' @export
plot.beset_zeroinfl <- function(object, metric = "MCE", se = TRUE){
  train <- dplyr::select(object$fit_stats, n_count_pred, n_zero_pred, form,
                         dplyr::starts_with(metric))
  names(train)[4] <- "train"
  xval <- dplyr::select(object$xval_stats, n_count_pred, n_zero_pred, form,
                        dplyr::starts_with(metric))
  names(xval)[4:5] <- c("cv", "cv_se")
  test <- try(dplyr::select(object$test_stats, n_count_pred, n_zero_pred, form,
                            dplyr::starts_with(metric)), silent = TRUE)
  if(class(test) == "try-error"){
    test <- NULL
  } else {
    names(test)[3] <- "test"
  }
  data <- dplyr::left_join(train, xval,
                           by = c("n_count_pred", "n_zero_pred", "form"))
  if(!is.null(test)) data <- dplyr::left_join(data, test)
  data$cv_lower <- with(data, cv - cv_se)
  data$cv_upper <- with(data, cv + cv_se)
  xmax <- max(data$n_count_pred)
  color_legend <- c("Train" = "grey", "CV" = "red", "Test" = "blue")
  y_lab <- switch(metric,
                  MCE = ylab("Mean Cross Entropy (-loglik / N)"),
                  MSE = ylab("Mean Squared Error"),
                  R2 = ylab(bquote(~R^2)))
  p <- ggplot(data = data) + facet_wrap(~n_zero_pred) + theme_bw() +
    ggtitle(title) + xlab("Number of Predictors of Count Component") + y_lab +
    scale_x_continuous(breaks = 0:xmax) +
    scale_color_manual(name = "", values = color_legend) +
    geom_line(aes(x = n_count_pred, y = train, color = "Train")) +
    geom_line(aes(x = n_count_pred, y = cv, color = "CV")) +
    ggtitle(paste("Numbers at top of each panel indicate",
                "\nnumber of predictors of zero-inflation component.")) +
    theme(plot.title = element_text(hjust=0, size=10))

  if(se){
    p <- p +
      geom_errorbar(aes(x = n_count_pred, ymin = cv_lower, ymax = cv_upper),
                    width = 0.2, color = "red")
  }
  if(!is.null(test)){
    p <- p + geom_line(aes(x = n_count_pred, y = test_CE, color = "Test"))
  }
  print(p)
}
