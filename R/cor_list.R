#' Construct a Correlation List
#'
#' Construct a \code{cor_list} object from a correlation matrix.
#'
#' The R stats package returns bivariate correlations as a matrix of values,
#' which can be difficult to parse when there are many correlations.
#' \code{cor_list} converts the correlation matrix into an adjacency list data
#' structure, removing the diagonal entries.
#'
#' @return A \code{cor_list} object with the following three vectors:
#'  \describe{
#'    \item{\code{i}}{the row names of the correlation matrix}
#'    \item{\code{j}}{the column names of the correlatin matrix}
#'    \item{\code{coef}}{the correlation coefficient corresponding to element
#'    \code{(i, j)} of the correlation matrix}
#'  }
#'
#' @param cor_matrix Correlation matrix
#'
#' @export

cor_list <- function(cor_matrix){
  i <- rep(rownames(cor_matrix), times = ncol(cor_matrix))
  j <- rep(colnames(cor_matrix), each = nrow(cor_matrix))
  coef = as.vector(cor_matrix)
  output <- data.frame(i = i, j = j, coef = coef, stringsAsFactors = FALSE)
  structure(dplyr::filter(output, i != j), class = "cor_list")
}
