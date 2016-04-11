#' Plot Method for the \code{beset_elnet} Class
#'
#' \code{plot.beset_elnet} takes the output of a \code{\link{beset_elnet}}
#'  object and creates a line plot using \code{\link[ggplot2]{ggplot}}.
#'
#' \code{plot.beset_elnet} produces a plot showing cross-validation error as a
#' function of alpha (L1-L2 mixing parameter) and lambda (regularization
#' parameter) from a cross-validated elastic net regression, as implemented by
#'  \code{\link{beset_elnet}}.
#'
#' @param object An object of class \code{beset_elnet}.
#'
#' @param title Optional character vector giving plot title.
#'
#' @import ggplot2
#'
#' @export

plot.beset_elnet <- function(object, title = ""){
  data <- object$results
  data$alpha <- as.factor(data$alpha)
  ggplot(data = data) +
    theme_bw() +
    ggtitle(title) +
    xlab("Regularization parameter") +
    ylab("Cross-validation error") +
    scale_x_log10() +
    geom_ribbon(aes(x = lambda, ymin = cve_lo,
                    ymax = cve_hi, color = alpha), fill = "grey") +
    labs(color = "Mixing parameter")
}

