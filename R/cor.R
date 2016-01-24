#' Correlation
#'
#' This function overrides the \code{\link[stats]{cor}} function provided by the
#' stats library to give the option of returning correlations as a
#' \code{\link{cor_list}} object.
#'
#' \code{cor} differs from R's native \code{\link[stats]{cor}} only in that it
#' uses pairwise-complete observations and returns a \code{\link{cor_list}}
#' object by default. If you want \code{beset::cor} to mimic the default
#' behavior of \code{stats::cor}, set \code{as_matrix = TRUE} and
#' \code{use = "everything"}.
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
#' @seealso \code{\link[stats]{cor}}, \code{\link{cor_list}}
#'
#' @export

cor <- function(x, y = NULL, as_matrix = FALSE, use = "pairwise.complete.obs",
                     method = "pearson"){
  output <- stats::cor(x, y = y, use = use, method = method)
  stat <- "r"
  if(method != "pearson") stat <- ifelse(method == "kendall", "tau", "rho")
  if(!as_matrix){
    output <- cor_list(output)
    attr(output, "row.names") <- NULL
  }
  attr(output, "coef") <- stat
  output
}

