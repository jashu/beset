#' Create Folds for Cross-Validation
#'
#' \code{create_folds} randomly splits the data into \code{k} folds, with
#' optional repetition of the splitting procedure under different randomization
#' seeds
#'
#' \code{create_folds} is used by \code{beset} modeling functions to assign
#' folds, stratified by the response variable (\code{y}). When \code{y} is a
#' factor, random fold assignment is done within the levels of \code{y} to
#' balance the class distribution across folds. When \code{y} is numeric, the
#' sample is split into strata based on percentiles and random fold assignment
#' is done within strata. The number of groups is set dynamically based on
#' sample size and \code{k}. In addition, if \code{y} is zero-inflated, the
#' algorithm attempts to balance the proportion of \code{0} values across folds.
#' Optionally, if the \code{reps} argument is set to a value greater than
#' \code{1}, the function repeats the stratified fold assignment by this number
#' of \code{reps}, using different random seeds each time.
#'
#' @return A list of \eqn{k * reps} row-position integers corresponding to the
#' training data for each independent fold of each repetition. (Test data are
#' implied by the rows that are omitted.) The names of the list elements give
#' the fold membership using the pattern "Fold\emph{i}.Rep\emph{j}".
#'
#' @section Note:
#' Previous versions of \code{beset} imported the function
#' \code{\link[caret]{createMultiFolds}} to set up repeated k-fold
#' cross-validation. While the object returned by \code{create_folds} copies the
#' structure of what would be returned by \code{\link[caret]{createMultiFolds}},
#' the algorithm for fold assignment is different. This may result in slightly
#' different cross-validation results if you attempt to reproduce an analysis
#' performed with a previous version.
#'
#' @param y A vector of outcomes.
#'
#' @param k An integer for the number of folds (default = 10).
#'
#' @param reps An integer for the number of repetitions (default = 1).
#'
#' @param seed An integer used to reproduce the random fold assignments.
#'
#' @seealso \code{\link{data_partition}}
#'
#' @export

create_folds <- function(y, k = 10, reps = 1, seed = 42){
  set.seed(seed)
  seed_list <- sample.int(k*reps, size = reps)
  folds <- purrr::map(seed_list, ~ stratify_folds(y, k, .x))
  fold_ids <- expand.grid(Fold = 1:k, Rep = 1:reps)
  fold_names <- paste("Fold", fold_ids$Fold, ".Rep", fold_ids$Rep, sep = "")
  fold_ids <- purrr::map2(fold_ids$Fold, fold_ids$Rep, ~ which(folds[[.y]] != .x))
  names(fold_ids) <- fold_names
  fold_ids
}


assign_folds <- function(y, k = 10, seed = 42){
  fold_ids <- rep(1:k, ceiling(length(y)/k))
  set.seed(seed)
  sample(fold_ids)[1:length(y)]
}

stratify_folds <- function(y, k = 10, seed = 42){
  folds <- vector("integer", length(y))
  if(is.numeric(y)){
    if(min(y) == 0){
      # make separate fold assignments for zero vs. non-zero values to insure
      # both 0 and count process is represented in every fold
      folds[y == 0] <- assign_folds(y[y == 0], k, seed)
      folds[y != 0] <- stratify_folds(y[y != 0], k, seed)
    } else {
      cuts <- floor(length(y)/k)
      if(cuts < 2) cuts <- 2
      if(cuts > 5) cuts <- 5
      breaks <- unique(quantile(y, probs = seq(0, 1, length = cuts)))
      y <- cut(y, breaks, include.lowest = TRUE) %>% as.integer
    }
  }
  purrr::walk(unique(y), function(x){
    folds[y == x] <<- assign_folds(y[y == x], k, seed)
  })
  folds
}


