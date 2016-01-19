#' Correlation List
#'
#' \code{cor_list} is a wrapper for \code{\link[stats]{cor}} that returns
#' correlations in a list (as opposed to a matrix)
#'
#' The R stats package returns bivariate correlations as a matrix of values,
#' which can be difficult to read when there are many correlations and which
#' contains redundant and useless information; i.e., one only needs to look at
#' values either above or below the diagonal. \code{cor_list} extracts unique
#' pairwise combinations of variables and displays their correlations in an
#' easy-to-read list. It also uses pairwise complete observations by default,
#' which is the more common use case for exploratory analyses, but is otherwise
#' equivalent to \code{\link[stats]{cor}} in behavior.
#'
#' @param data A numeric matrix or data frame containing only numeric variables.
#'
#' @param use	An optional character string giving a method for computing
#' covariances in the presence of missing values. This must be (an abbreviation
#' of) one of the strings "everything", "all.obs", "complete.obs",
#' "na.or.complete", or "pairwise.complete.obs".
#'
#' @param  method	A character string indicating which correlation coefficient
#' is to be computed. One of "pearson" (default), "kendall", or "spearman": can
#' be abbreviated.
#'
#' @return A numeric vector of correlation coefficients, named according to each
#'  unique pairwise combination of variables.
#'
#' @seealso \code{\link{cor}}
#'
#' @export

cor_list <- function(data, use = "pairwise.complete.obs",
                     method = "pearson"){
  cor_matrix <- cor(data, use = use, method = method)
  upper_tri <- upper.tri(cor_matrix)
  r <- as.vector(cor_matrix)[as.vector(upper_tri)]
  indices <- which(upper_tri, arr.ind = TRUE)
  from <- rownames(cor_matrix)[indices[,1]]
  to <- colnames(cor_matrix)[indices[,2]]
  names(r) <- paste(from, to, sep = "<->")
  stat <- "r"
  if(method != "pearson"){
    stat <- ifelse(method == "kendall", "tau", "rho")
  }
  attr(r, "statistic") <- stat
  return(r)
}

