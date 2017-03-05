#' Beset GLM with Elasticnet Regularization
#'
#' \code{beset_elnet} is a wrapper to \code{\link[glmnet]{glmnet}} for fitting
#'  generalized linear models via penalized maximum likelihood, providing
#'  automated data preprocessing and selection of both the elasticnet penalty
#'  and regularization parameter through repeated cross-validation.
#'
#' @param form A model \code{\link[stats]{formula}}.
#'
#' @param train_data Data frame with the variables in \code{form} and the data
#' to be used for model fitting.
#'
#' @param test_data Optional data frame with the variables in \code{form} and
#' the data to be used for model validation.
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
#' @param standard_coef Logical flag to return standardized regression
#' coefficients. Default is \code{standard_coef = TRUE}. Note that this refers
#' to how the coefficients are returned; the elastic net requires all \code{x}
#' variables to be at the same scale, so if \code{standard_coef = FALSE}
#' these variables will still be standardized for purposes of model fitting,
#' but the coefficients will be converted back to original scale before they are
#' returned.
#'
#' @param ... Additional parameters to be passed to
#' \code{\link[glmnet]{glmnet}}.
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

beset_elnet <- function(form, train_data, test_data = NULL,
                        family = "gaussian", alpha = c(.05, .5, .95),
                        standard_coef = TRUE, ...,
                        n_folds = 10, n_repeats = 10,
                        seed = 42, n_cores = 2){
  #==================================================================
  # Create model frame and extract response name and vector
  #------------------------------------------------------------------
  mf <- stats::model.frame(form, data = train_data, na.action = stats::na.omit)
  if(!is.null(test_data)){
    test_data <- model.frame(form, data = test_data, na.action = stats::na.omit)
  }
  n_drop <- nrow(train_data) - nrow(mf)
  if(n_drop > 0)
    warning(paste("Dropping", n_drop, "rows with missing data."),
            immediate. = TRUE)
  if(standard_coef) mf <- dplyr::mutate_if(mf, is.numeric, scale)
  y <- mf[,1]
  if(grepl("binomial", family)) y <- as.factor(y)
  x <- stats::model.matrix(form, data = mf)[,-1]

  #======================================================================
  # Obtain fit statistics
  #----------------------------------------------------------------------
  fits <- lapply(alpha, function(a){
    glmnet::glmnet(x, y, family, alpha = a)#, ...)
   })
  names(fits) <- alpha
  lambda_seq <- sapply(fits, function(x) unique(round(x$lambda, 3)))
  lambda_length <- sapply(lambda_seq, function(x) length(x))
  alpha_seq <- mapply(function(a, l) rep(a, l), a = alpha, l = lambda_length)
  fit_stats <- data.frame(alpha = unlist(alpha_seq),
                          lambda = unlist(lambda_seq))
  metrics <- mapply(function(fit, lambda){
    y_hat <- stats::predict(fit, x, lambda, type = "response")
    apply(y_hat, 2, function(x) predict_metrics_(y, x, family))
  }, fit = fits, lambda = lambda_seq)
  fit_stats$MCE <- unlist(sapply(metrics, function(x)
    sapply(x, function(s) s$mean_cross_entropy)))
  fit_stats$MSE <- unlist(sapply(metrics, function(x)
    sapply(x, function(s) s$mean_squared_error)))
  fit_stats$R2 <- unlist(sapply(metrics, function(x)
    sapply(x, function(s) round(s$R_squared, 3))))

  #==================================================================
  # Obtain independent test stats
  #------------------------------------------------------------------
  test_stats <- NULL
  if(!is.null(test_data)){
    test_stats <- data.frame(alpha = unlist(alpha_seq),
                             lambda = unlist(lambda_seq))
    if(standard_coef){
      test_data <- dplyr::mutate_if(test_data, is.numeric, scale)
    }
    y_test <- test_data[,1]
    if(grepl("binomial", family)) y_test <- factor(y_test)
    x_test <- model.matrix(form, data = test_data)[,-1]
    metrics <- mapply(function(fit, lambda){
      y_hats <- predict(fit, x_test, lambda, type = "response")
      apply(y_hats, 2, function(y_hat) predict_metrics_(y_test, y_hat, family))
    }, fit = fits, lambda = lambda_seq)
    test_stats$MCE <- unlist(sapply(metrics, function(x)
      sapply(x, function(s) s$mean_cross_entropy)))
    test_stats$MSE <- unlist(sapply(metrics, function(x)
      sapply(x, function(s) s$mean_squared_error)))
    test_stats$R2 <- unlist(sapply(metrics, function(x)
      sapply(x, function(s) round(s$R_squared, 3))))
  }

  #==================================================================
  # Obtain cross-validation stats
  #------------------------------------------------------------------
  set.seed(seed)
  fold_ids <- caret::createMultiFolds(y, k = n_folds, times = n_repeats)
  fold_list <- rep(fold_ids, length(alpha))
  alpha_list <- rep(alpha, each = n_folds * n_repeats)
  lambda_list <- rep(lambda_seq, each = n_folds * n_repeats)
  cl <- parallel::makeCluster(n_cores)
  parallel::clusterExport(cl, "predict_metrics_")
  metrics <- parallel::clusterMap(
    cl, function(i, alpha, lambda, predictors, response, ...){
      fit <- glmnet::glmnet(x = predictors[i,], y = response[i], alpha = alpha,
                            ...)
      y <- response[-i]
      y_hat <- predict(fit, predictors[-i,], lambda, "response")
      apply(y_hat, 2, function(x) predict_metrics_(y, x, "gaussian"))
      }, i = fold_list, alpha = alpha_list, lambda = lambda_list,
    MoreArgs = list(predictors = x, response = y, family = family)#, ...)
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
                                 MCE_SE = get_se(MCE),
                                 MSE_SE = get_se(MSE),
                                 R2_SE = get_se(R2),
                                 MCE = mean(MCE),
                                 MSE = mean(MSE),
                                 R2 = mean(R2))
  cv_stats <- dplyr::ungroup(cv_stats)
  cv_stats <- dplyr::select(cv_stats, alpha, lambda,
                              MCE, MCE_SE, MSE, MSE_SE, R2, R2_SE)
  cv_stats <- dplyr::arrange(cv_stats, alpha, desc(lambda))

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

rec_elnet <- function(x, y, ..., n_vars = NULL, threshold = NULL, seed = 42){
  n_observed <- nrow(x)
  n_features <- ncol(x)
  snp_name <- colnames(x)
  if(is.null(n_vars)) n_vars <- ceiling(sqrt(n_features))
  if(is.null(threshold)) threshold <- 1 / n_observed
  n_samples <- n_vars * 10
  betas <- matrix(NA, nrow = n_samples, ncol = n_features)
  intercepts <- vector(mode = "numeric", length = n_samples)
  set.seed(seed, kind = "default")
  pb <- txtProgressBar(min = 0, max = n_samples, style = 3)
  for(i in 1:n_samples){
    importance <- apply(betas, 2, mean, na.rm = T)
    importance[is.nan(importance)] <- 1
    importance <- abs(importance)
    importance[importance < threshold] <- 0
    rows_in_train <- sample.int(nrow(x), replace = TRUE)
    cols_in_train <- sample.int(n_features, n_vars, prob = importance)
    temp_x <- x[rows_in_train, cols_in_train]
    temp_y <- y[rows_in_train]
    temp_pf <- 1 / importance[cols_in_train]
    temp_mod <- gcdnet::cv.gcdnet(temp_x, temp_y, lambda2 = 1,
                                  pf = temp_pf, ...)
    betas[i, cols_in_train] <- as(coef(temp_mod), "matrix")[-1,]
    intercepts[i] <- as(coef(temp_mod), "matrix")[1,]
    setTxtProgressBar(pb, i)
  }
  importance <- apply(betas, 2, mean, na.rm = T)
  importance[is.nan(importance)] <- 0
  importance <- abs(importance)
  importance[importance < threshold] <- 0
  structure(
    list(
      cv_gcdnet = gcdnet::cv.gcdnet(x[, importance > 0], y, lambda2 = 1,
                                    pf = 1 / importance[importance > 0], ...),
      b0_bag = intercepts,
      beta_bag = as(betas, "dgCMatrix")),
    class = "rr_elnet"
  )
}

