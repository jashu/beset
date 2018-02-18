#' Beset GLM with Elasticnet Regularization
#'
#' \code{beset_elnet} is a wrapper to \code{\link[glmnet]{glmnet}} for fitting
#'  generalized linear models via penalized maximum likelihood, providing
#'  automated data preprocessing and selection of both the elasticnet penalty
#'  and regularization parameter through repeated k-fold cross-validation.
#'
#' @param form A model \code{\link[stats]{formula}}.
#'
#' @param data Either a \code{\link{data_partition}} object containing data sets
#' to be used for both model training and testing, or a single data frame that
#' will be used for model training only.
#'
#' @param test_data Optional data frame with the variables in \code{form} and
#' the data to be used for model validation. Not necessary and ignored if a
#' \code{\link{data_partition}} object is passed via the \code{data} parameter.
#'
#' @param family Character string naming the error distribution to be used in
#' the model. Currently supported options are \code{"gaussian"} (default),
#' \code{"binomial"}, and \code{"poisson"}.
#'
#' @param alpha Numeric vector of alpha values between 0 and 1 to use as tuning
#' parameters. \code{alpha = 0} results in ridge regression, and \code{alpha =
#' 1} results in lasso regression. Values in between result in a mixture of L1
#' and L2 penalties. (Values closer to 0 weight the L2 penalty more heavily,
#' and values closer to 1 weight the L1 penalty more heavily.)
#'
#' @param n_folds Integer indicating the number of cross-validation folds.
#'
#' @param n_repeats Number of times cross-validation should be repeated.
#'
#' @param n_cores Number of cores to use for parallel execution. If not
#' specified, the number of cores is set to 2. To determine the theoretical
#' maximum number of cores you have available, see
#' \code{\link[parallel]{detectCores}}, but note that the actual number of cores
#' available may be less. See \code{\link[parallel]{parallel}} for
#' more information.
#'
#' @param seed Seed for random number generator, used for assigning observations
#' to cross-validation folds.
#'
#' @seealso \code{\link[glmnet]{glmnet}}
#'
#' @export

beset_elnet <- function(form, data, test_data = NULL,
                        family = "gaussian", alpha = c(.01, .5, .99), ...,
                        n_folds = 10, n_repeats = 10,
                        seed = 42, n_cores = 2){
  #==================================================================
  # Check if data is data_partition object and set up accordingly
  #------------------------------------------------------------------
  if(inherits(data, "data_partition")){
    test_data <- data$test
    data <- data$train
  }

  #==================================================================
  # Create model frame and extract response name and vector
  #------------------------------------------------------------------
  mf <- model_frame(form, data)
  # Correct non-standard column names
  names(mf) <- make.names(names(mf))
  # Do the same for test_data if it exists
  if(!is.null(test_data)){
    test_data <- model_frame(form, data = test_data)
    names(test_data) <- make.names(names(test_data))
    if(!all(names(mf) %in% names(test_data)))
      stop("'test_data' must contain same variables as 'data'")
  }
  # Warn user if any rows were dropped
  n_drop <- nrow(data) - nrow(mf)
  if(n_drop > 0){
    warning(paste("Dropping", n_drop, "rows with missing data.\n"),
            immediate. = TRUE)
  }
  # cache name of the response variable
  response <- names(mf)[1]
  # extract response vector
  y <- mf[[1]]
  # insure y is a factor if family is binomial
  if(family == "binomial" && !is.factor(y)) y <- factor(y)

  #======================================================================
  # Obtain fit statistics
  #----------------------------------------------------------------------
  x <- stats::model.matrix(form, data = mf)[, -1]
  if(!is.null(test_data)){
    x_test <- stats::model.matrix(form, data = test_data)[, -1]
    x_test <- try(x_test[, colnames(x)], silent = TRUE)
    if(inherits(x_test, "try-error")){
      stop("Train and test data have different predictors or predictor levels")
    }
    y_test <- test_data[[1]]
  }
  fits <- lapply(alpha, function(a){
    glmnet::glmnet(x, y, family, alpha = a)
   })
  names(fits) <- alpha
  lambda_seq <- sapply(fits, function(x) unique(round(x$lambda, 3)))
  lambda_length <- sapply(lambda_seq, function(x) length(x))
  alpha_seq <- mapply(function(a, l) rep(a, l), a = alpha, l = lambda_length)
  fit_stats <- data.frame(alpha = unlist(alpha_seq),
                          lambda = unlist(lambda_seq))
  metrics <- mapply(function(fit, lambda){
    y_hats <- stats::predict(fit, x, lambda, type = "response")
    y_obs <- if(is.factor(y)) as.integer(y) - 1 else y
    apply(y_hats, 2, function(y_hat) predict_metrics_(y_obs, y_hat, family))
  }, fit = fits, lambda = lambda_seq)
  fit_stats$mce <- unlist(sapply(metrics, function(x)
    sapply(x, function(s) s$mean_cross_entropy)))
  fit_stats$mse <- unlist(sapply(metrics, function(x)
    sapply(x, function(s) s$mean_squared_error)))
  fit_stats$r2 <- unlist(sapply(metrics, function(x)
    sapply(x, function(s) round(s$R_squared, 3))))

  #==================================================================
  # Obtain independent test stats
  #------------------------------------------------------------------
  test_stats <- NULL
  if(!is.null(test_data)){
    test_stats <- data.frame(alpha = unlist(alpha_seq),
                             lambda = unlist(lambda_seq))
    metrics <- mapply(function(fit, lambda){
      y_hats <- stats::predict(fit, x_test, lambda, type = "response")
      y_obs <- if(is.factor(y_test)) as.integer(y_test) - 1 else y_test
      apply(y_hats, 2, function(y_hat) predict_metrics_(y_obs, y_hat, family))
    }, fit = fits, lambda = lambda_seq)
    test_stats$mce <- unlist(sapply(metrics, function(x)
      sapply(x, function(s) s$mean_cross_entropy)))
    test_stats$mse <- unlist(sapply(metrics, function(x)
      sapply(x, function(s) s$mean_squared_error)))
    test_stats$r2 <- unlist(sapply(metrics, function(x)
      sapply(x, function(s) round(s$R_squared, 3))))
  }

  #==================================================================
  # Obtain cross-validation stats
  #------------------------------------------------------------------
  fold_ids <- create_folds(y, n_folds, n_repeats, seed)
  fold_list <- rep(fold_ids, length(alpha))
  alpha_list <- rep(alpha, each = n_folds * n_repeats)
  lambda_list <- rep(lambda_seq, each = n_folds * n_repeats)
  cl <- parallel::makeCluster(n_cores)
  lib_paths <- .libPaths()
  parallel::clusterExport(cl, "lib_paths", envir = environment())
  parallel::clusterEvalQ(cl, .libPaths(lib_paths))
  metrics <- parallel::clusterMap(
    cl, function(i, alpha, lambda, x, y, family){
      fit <- glmnet::glmnet(x = x[i,], y = y[i], alpha = alpha,
                            family = family)
      y_hat <- stats::predict(fit, x[-i,], lambda, "response")
      y <- if(is.factor(y)) as.integer(y) - 1
      apply(y_hat, 2, function(x) predict_metrics_(y[-i], x, family))
      }, i = fold_list, alpha = alpha_list, lambda = lambda_list,
    MoreArgs = list(x = x, y = y, family = family)
  )
  parallel::stopCluster(cl)
  cv_stats <- data.frame(
    fold_id = rep(names(fold_ids), each = nrow(fit_stats)),
    alpha = rep(fit_stats$alpha, n_repeats*n_folds),
    lambda = rep(fit_stats$lambda, n_repeats*n_folds))
  cv_stats <- dplyr::arrange(cv_stats, alpha, fold_id, dplyr::desc(lambda))
  cv_stats$MCE <- unlist(sapply(metrics, function(lambdas)
    sapply(lambdas, function(x) x$mean_cross_entropy)))
  cv_stats$MSE <- unlist(sapply(metrics, function(lambdas)
    sapply(lambdas, function(x) x$mean_squared_error)))
  cv_stats$R2 <- unlist(sapply(metrics, function(lambdas)
    sapply(lambdas, function(x) round(x$R_squared, 3))))
  get_se <- function(x){
    sd(x)/sqrt(length(x))
  }
  cv_stats <- dplyr::group_by(cv_stats, alpha, lambda)
  cv_stats <- dplyr::summarize(cv_stats,
                                 mce_SE = get_se(MCE),
                                 mse_SE = get_se(MSE),
                                 r2_SE = get_se(R2),
                                 mce = mean(MCE),
                                 mse = mean(MSE),
                                 r2 = mean(R2))
  cv_stats <- dplyr::ungroup(cv_stats)
  cv_stats <- dplyr::select(cv_stats, alpha, lambda,
                              mce, mce_SE, mse, mse_SE, r2, r2_SE)
  cv_stats <- dplyr::arrange(cv_stats, alpha, dplyr::desc(lambda))

  #======================================================================
  # Construct beset_glm object
  #----------------------------------------------------------------------
  structure(list(stats = list(fit = fit_stats, cv = cv_stats,
                              test = test_stats),
                 cv_params = list(n_folds = n_folds, n_repeats = n_repeats,
                                    seed = seed, fold_ids = fold_ids),
                 model_fits = fits,
                 model_params = list(x = x, y = y, family = family)),
            class = "beset_elnet")

}

