#' Leave-One-Out Correlations
#'
#' \code{cor_loo} computes the leave-one-out estimate for all Pearson
#' correlations among the columns of a matrix or data frame.
#'
#' The leave-one-out (LOO) cross-validation correlation provides an estimate of
#' Pearson's \eqn{r} that corrects for inflation in the absolute magnitutde of
#' \eqn{r} associated with small sample sizes, in which single obvervations can
#' have undue influence on the size of the correlation coefficient. It is
#' equivalent to fitting a linear model to all possible sets of \eqn{n-1}
#' observations, predicting the left-out observation, and using the root mean
#' square error (RMSE) of the predicitions to recalculate \eqn{r}. With linear
#' models, this quantity can be computed without repeating the leave-one-out
#' procedure \eqn{N} times by using the companion \code{\link{predict_R2}}
#' function, which computes the predicted \eqn{r^2}{r-squared} using the PRESS
#' statistic. The LOO estimate of \eqn{r} is then given as the square root of
#' the predicted \eqn{r^2}{r-squared} for each bivariate linear regression
#' model, with the sign of the LOO \eqn{r} assigned to match that of the
#' corresponding \eqn{beta} coefficient.
#'
#' @section Note:
#' Correlations are computed using pairwise complete observations.
#'
#' @param data A numeric matrix or data frame containing only numeric variables.
#'
#' @return A numeric vector of correlation coefficients, named according to each
#'  unique pairwise combination of variables.
#'
#' @seealso \code{\link{predict_R2}}
#'
#' @export

cor_loo <- function(data){
  cor_matrix <- cor(data, use = "pairwise")
  upper_tri <- upper.tri(cor_matrix)
  indices <- which(upper_tri, arr.ind = TRUE)
  from <- rownames(cor_matrix)[indices[,1]]
  to <- colnames(cor_matrix)[indices[,2]]
  indices <- lapply(seq_len(nrow(indices)), function(i) indices[i,])
  loo_cor <- function(i){
    x <- unlist(data[,i[1]])
    y <- unlist(data[,i[2]])
    model <- lm(y ~ x)
    pR2 <- predict_R2(model)
    if(!is.na(pR2) && pR2 < 0) pR2 <- 0
    loo_r <- ifelse(coef(model)[2] < 0, sqrt(pR2) * -1, sqrt(pR2))
  }
  loo_r <- parallel::mcmapply(loo_cor, indices)
  names(loo_r) <- paste(from, to, sep = "<->")
  return(loo_r)
}
