#' Cross-Validation Mean
#'
#' Estimate cross-validation performance that would be obtained by averaging
#' predictions of more than one model.
#'
#' @param x A "cross_valid" object returned by \code{\link{validate}}
#'
#' @param ... Additional "cross_valid" objects.
#'
#' @export

mean.cross_valid <- function(x, ...){
  ob <- list(x, ...)
  fold_ids <- map(ob, "fold_assignments")
  ids <- fold_ids[[1]]
  if(!all(map_lgl(fold_ids, ~identical(ids, .x)))){
    stop("Models must be cross-validated using identical fold assignments.")
  }
  y <- as.numeric(ob[[1]]$parameters$y)
  family <- ob[[1]]$parameters$family
  n_folds <- ob[[1]]$parameters$n_folds
  fold_ids <- map(ids, function(x) map(1:n_folds, ~ which(x == .x)))
  predictions <- map(ob, ~ as.matrix(.x$predictions)) %>% reduce(`+`)
  predictions <- as.data.frame(predictions / length(ob))
  stats_per_fold <- map2(
    predictions, fold_ids, function(yhat, i){
      map(i, ~ predict_metrics_(y[.x], yhat[.x], family))
    }
  ) %>% map(~transpose(.x) %>% simplify_all) %>% transpose %>% simplify_all
  stats_per_rep <- map(predictions, ~ predict_metrics_(y, .x, family)) %>%
    transpose %>% simplify_all
  ensemble_stats <- list(
    mean = map(stats_per_rep, mean),
    btwn_fold_se = map(stats_per_fold, ~ sd(.x) / sqrt(10)),
    btwn_rep_range = map(stats_per_rep, range)
  ) %>% transpose
  structure(
    list(
      stats = ensemble_stats,
      predictions = predictions,
      fold_assignments = ids,
      parameters = ob[[1]]$parameters
    ),
    class = class(ob[[1]])
  )
}
