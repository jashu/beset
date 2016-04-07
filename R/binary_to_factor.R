#' Convert Binary Logical or Numeric Variable to Factor
#'
#' Converts binary logical or numeric vectors to factor objects, assigning the
#' positive instance (\code{TRUE} or \code{1}) to be the reference level.
#'
#' Functions in the
#' \code{\href{http://topepo.github.io/caret/index.html}{caret}}
#' package expect binary variables to be encoded as factors, with the positive
#' instance of the class assigned to the first (reference) level of the factor.
#' If your class labels are coded as numeric binary (0s and 1s), a simple
#' conversion to factor without reordering the levels will reverse the labels
#' because factor levels are assigned based on alphanumeric order. Thus, 0 will
#' be treated as the positive class and 1 as the negative class. The same is
#' true for class labels coded as logical binary (\code{FALSE} and \code{TRUE}).
#' This does not affect the model fit itself, but it will cause classification
#' statistics to be reversed, e.g., sensitivity will be labeled specificity and
#' vice versa.
#'
#' \code{binary_to_factor} is a wrapper for \code{\link[base]{factor}} with
#' automatic reordering of factor levels and a simpler method of specifying
#' factor labels.
#'
#' @param binary Logical vector or numeric vector containing only 0s and 1s.
#'
#' @param neg Optional character string of length 1 with a label for the
#' negative class (the class coded as \code{FALSE} or \code{0}).
#'
#' @param pos Optional character string of length 1 with a label for the
#' positive class (the class coded as \code{TRUE} or \code{1}).
#'
#' @return an object of class "\code{factor}"
#'
#' @export

binary_to_factor <- function(binary, neg = "class_0", pos = "class_1"){
  binary <- as.numeric(binary)
  factor(binary, levels = c(1,0), labels = c(pos, neg))
}
