#' Correlation
#'
#' This function overrides the \code{\link[stats]{cor}} function provided by the
#' stats library to give the option of returning correlations as a
#' \code{\link{cor_list}} object.
#'
#' \code{cor} differs from R's native \code{\link[stats]{cor}} only in that it
#' uses pairwise-complete observations and returns a \code{\link{cor_list}}
#' object by default, which has special summary methods that can be used to
#' selectively view the relationships of interest. See examples and
#' \code{\link{summary.cor_list}}.
#'
#' @section Note:
#' If you want \code{beset::cor} to mimic the default
#' behavior of \code{stats::cor}, set \code{as_matrix = TRUE} and
#' \code{use = "everything"}.
#'
#' @examples
#' # Create a correlation list for the numeric variables from the iris data set
#' iris_cors <- cor(iris[,-5])
#'
#' # Look at all correlations with Sepal.Length
#' summary(iris_cors, x1 = Sepal.Length)
#'
#' # Look at all correlations between sepal measurements and petal measurements
#' summary(iris_cors, x1 = starts_with("Sepal"), x2 = starts_with("Petal"))
#'
#' # Look at all correlations with Sepal.Length, excluding Sepal.Width
#' summary(iris_cors, x1 = Sepal.Length, x2 = -Sepal.Width)
#'
#' # Look at all correlations in their original order
#' summary(iris_cors, sort = FALSE)
#'
#' @inheritParams stats::cor
#'
#' @param as_matrix Logical indicating whether correlations should be returned
#' as a matrix.
#'
#' @return If \code{as_matrix = FALSE} (the default), a \code{\link{cor_list}}
#' object. If \code{as_matrix = TRUE}, a correlation matrix. Both are assigned
#' the attribute "coef" which indicates the statistic that is being returned,
#' e.g., Pearson's \emph{r}, Kendall's \emph{tau}, or Spearman's \emph{rho}.
#'
#' @seealso \code{\link[stats]{cor}}, \code{\link{cor_list}},
#' \code{\link{summary.cor_list}}
#'
#' @export

cor <- function(x, y = NULL, as_matrix = FALSE, use = "pairwise.complete.obs",
                     method = "pearson"){
  output <- stats::cor(x, y = y, use = use, method = method)
  if (class(output) == "numeric"){
    x_name <- names(x)[1]
    if (is.null(x_name)) x_name <- "x"
    y_name <- names(y)
    if (is.null(y_name)) y_name <- names(x)[2]
    if (is.null(y_name)) y_name <- "y"
    output <- structure(list(
      x1 = c(x_name, y_name),
      x2 = c(y_name, x_name),
      coef = rep(output, 2)),
      class = "cor_list")
  } else if (!as_matrix){
      output <- cor_list(output)
      attr(output, "row.names") <- NULL
  }
  stat <- "r"
  if(method != "pearson") stat <- ifelse(method == "kendall", "tau", "rho")
  attr(output, "coef") <- stat
  output
}

