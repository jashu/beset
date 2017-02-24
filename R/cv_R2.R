#' Cross-Validated R-Squared
#'
#' Returns a cross-validated estimate of R-squared, or how well the fitted model
#' object is expected to explain new data.
#'
#' @param n_cores Integer value indicating the number of workers to run in
#' parallel during subset search and cross-validation. By default, this will
#' be set to 2. You may wish to increase this depending on your hardware and OS.
#' See \code{\link[parallel]{parallel-package}} for more information.
#'
#' @param n_folds Integer indicating the number of folds to use for
#' cross-validation.
#'
#' @param n_repeats Integer indicating the number of times cross-validation
#' should be repeated.
#'
#' @param seed An integer used to seed the random number generator when
#' assigning observations to folds.
#'
#' @export
cv_r2 <- function(object, n_cores = 2, n_folds = 10, n_repeats = 10, seed = 42){
  model_type <- class(object)[1]
  y <- object$model[,1]
  set.seed(seed)
  fold_ids <- caret::createMultiFolds(y, k = n_folds, times = n_repeats)
  cl <- parallel::makeCluster(n_cores)
  r2 <- parallel::parSapply(cl, fold_ids, function(x, object){
    model <- switch(
      model_type,
      lm = lm(object$terms, data = object$model[x,]),
      glm = glm(object$terms, data = object$model[x,],
                family = object$family$family),
      negbin = MASS::glm.nb(object$terms, data = object$model[x,]),
      zeroinfl = pscl::zeroinfl(object$terms$full,
                                data = object$model[x,],
                                dist = object$dist)
    )
    beset::prediction_metrics(model, object$model[-x,])$R_squared
  }, object = object)
  boot_r2 <- boot::boot(r2, function(x, i) median(x[i], na.rm = TRUE), 1000,
                        parallel = "snow", cl = cl)
  parallel::stopCluster(cl)
  boot_CI <- boot::boot.ci(boot_r2, type = "bca")
  structure(list(cv_R2 = median(r2, na.rm = TRUE), `95% CI` = boot_CI$bca[4:5],
            R2 = r2), class = "cv_R2")
}

