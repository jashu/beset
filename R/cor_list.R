#' Construct a Correlation List
#'
#' Construct a \code{cor_list} object from a correlation matrix.
#'
#' The R stats package returns bivariate correlations as a matrix of values,
#' which can be difficult to parse when there are many correlations.
#' \code{cor_list} converts the correlation matrix into a list of all
#' row-column pairs, removing the diagonal entries.
#'
#' @return A \code{cor_list} object with the following three vectors:
#'  \describe{
#'    \item{\code{x1}}{variables on one side of the correlation}
#'    \item{\code{x2}}{variables on the other side of the correlation}
#'    \item{\code{coef}}{the correlation coefficient between \code{x1} and
#'      \code{x2}}
#'  }
#'
#' @param cor_matrix Correlation matrix
#'
#' @export

cor_list <- function(cor_matrix){
  i <- rep(rownames(cor_matrix), times = ncol(cor_matrix))
  j <- rep(colnames(cor_matrix), each = nrow(cor_matrix))
  coef = as.vector(cor_matrix)
  output <- data.frame(x1 = j, x2 = i, coef = coef, stringsAsFactors = FALSE)
  structure(dplyr::filter(output, x1 != x2), class = "cor_list")
}
