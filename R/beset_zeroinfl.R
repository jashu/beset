#' Best Subset Selection for Zero-inflated Count Data Regression
#'
#' \code{beset_zeroinfl} performs best subset selection using repeated
#' cross-validation to find the optimal number of predictors for zero-inflated
#' regression models.
#'
#' \code{beset_zeroinfl} performs best subset selection for zero-inflated
#' count models, fitting both a count model (e.g., negative binomial regression)
#' and zero-inflation model (e.g., logistic regression). Given that there are
#' two models involved, the search for the best combination of predictors is
#' not completely exhaustive, which would require fitting \eqn{(2^p)^2} models,
#' e.g., just 10 predictors have over 1 million possible combinations. Instead,
#' best subset selection is first performed separately for each of the two
#' models, using only an intercept term for the other model. This narrows the
#' search to \eqn{2(2^p) + p^2} models (e.g., 2148 for \eqn{p = 10}): the
#' best models for each possible number of predictors of the count process
#' (using a constant to predict the zero process) crossed with the best models
#' for each possible number of predictors of the zero process (using a constant
#' to predict the count process). Thus, for every \eqn{m} in \eqn{0, 1, ... p}
#' best predictors of the count process and every \eqn{n} in \eqn{0, 1, ... p}
#' best predictors of the zero-inflation process, \code{beset_zeroinfl} uses
#' \eqn{k}-fold cross-validation to select the model with the \eqn{m * n} number
#' of predictors that minimizes prediction error, i.e., how well the best models
#' trained using \eqn{k - 1} folds predict the out-of-fold samples.
#'
#' @section Cross-validation details:
#' When randomly assigning observations to cross-validation folds, responses
#' with zero values are first allocated to equalize (as much as possible) the
#' incidence of zero observations across folds. The remaining non-zero values
#' are then randomly allocated within subgroups based on percentiles to insure
#' that the folds will also be balanced in terms of the count distribution.
#'
#' \code{beset_zeroinfl} enforces the "one in ten rule" for each component of
#' the model with respect to the sample size expected for \eqn{k - 1} folds,
#' i.e., there must be a minimum of 10 observations of zero across \eqn{k - 1}
#' folds for every predictor in the zero-inflation model, and there must be
#' a minimum of 10 observations of non-zero counts across \eqn{k - 1} folds for
#' every predictor in the count model. \code{beset_zeroinfl} will screen for
#' this and, if possible, alter your settings for \code{p_zero_max} and/or
#' \code{p_count_max} to insure an adequate training sample size for all model
#' fits. If this happens, a warning message will be issued informing you of
#' the new settings.
#'
#' @section Warnings:
#' \enumerate{
#'  \item \code{beset_zeroinfl} handles missing data by performing listwise
#'    deletion. No other options for handling missing data are provided at this
#'    time. The user is encouraged to deal with missing values prior to running
#'    this function. A warning message will be issued indicating how much
#'    data are being deleted.
#'  \item \code{beset_zeroinfl} is intended for use with additive models only.
#'    If you include interaction terms in the model formula passed to
#'    \code{form}, they will be ignored (e.g. \code{A * B * C} will be treated
#'    as \code{A + B + C} and \code{A + B + A:B} will be treated as
#'    \code{A + B}). An exhaustive search over the space of possible
#'    interactions and/or non-linear effects is computationally prohibitive, but
#'    I hope to offer a greedy search option in the future. In the meantime and
#'    in general, I would recommend the
#'    \href{https://cran.r-project.org/web/packages/earth/index.html}{MARS}
#'    technique.
#'  \item \code{beset_zeroinfl} can be very slow and memory intensive.
#'  Attempting to run with more than 10 variables in the model data frame is not
#'  recommended.
#' }
#'
#' @seealso \code{\link{beset_glm}},
#' \code{\link[pscl]{zeroinfl}}, \code{\link[base]{set.seed}},
#' \code{\link{predict_metrics}}
#'
#' @param form A model \code{\link[stats]{formula}}. Note that you cannot
#' specify a separate model formula for count and zero components as when
#' calling \code{\link[pscl]{zeroinfl}} directly. Instead use a single formula
#' that would equal or nest the maximum complexity desired for both components,
#' and use the \code{p_count_max} and \code{p_zero_max} arguments if you want
#' to impose differential constraints on the complexity of each component.
#'
#' @param data Data frame with the variables in \code{form} and the data
#' to be used for model fitting.
#'
#' @param test_data Optional data frame with the variables in \code{form} and
#' the data to be used for model validation.
#'
#' @param family Character string naming the count model family. Options are
#' \code{"poisson"} (default) or \code{"negbin"}. (Note a log link is always
#' used with the count model.)
#'
#' @param link Character string naming the link function in the binary
#' zero-inflation model. Options are \code{"logit"} (default), \code{"probit"},
#' \code{"cloglog"}, \code{"cauchit"}, or \code{"log"}. (Note a binomial family
#' is always used with the zero-inflation model).
#'
#' @param p_count_max Maximum number of predictors allowed in count model.
#'
#' @param p_zero_max Maximum number of predictors allowed in zero model.
#'
#' @param n_folds Integer indicating the number of folds to use for
#' cross-validation.
#'
#' @param n_repeats Integer indicating the number of times cross-validation
#' should be repeated.
#'
#' @param n_cores Integer value indicating the number of workers to run in
#' parallel during subset search and cross-validation. By default, this will
#' be set to 2. To determine the theoretical maximum number of cores you have
#' available, see \code{\link[parallel]{detectCores}}, but note that the actual
#' number of cores available may be less. See
#' \code{\link[parallel]{parallel-package}} for more information.
#'
#' @param seed An integer used to seed the random number generator when
#' assigning observations to folds.
#'
#' @importFrom utils combn
#'
#' @return A "beset_zeroinfl" object with the following components:
#' \enumerate{
#'  \item\describe{
#'    \item{params}{list of parameters used in function call that will be
#'    needed to reproduce results}
#'  }
#'   \item\describe{
#'     \item{model_data}{data frame extracted from \code{data} and used to
#'      identify best subsets.}
#'  }
#'  \item\describe{
#'    \item{stats}{A list of five data frames, each containing metrics
#'    describing model fits or predictions:
#'    \describe{\item{count_fit}{a data frame containing fit statistics for every
#'      possible combination of predictors of count data:
#'      \describe{
#'      \item{n_count_pred}{the number of count predictors in model; note that
#'       the number of predictors for a factor variable corresponds to the
#'       number of factor levels minus 1}
#'      \item{count_pred}{rhs of formula for count model}
#'      \item{aic}{\eqn{-2*log-likelihood + k*npar}, where \eqn{npar} represents
#'      the number of parameters in the fitted model, and \eqn{k = 2}}
#'      \item{dev}{twice the difference between the log-likelihoods of the
#'              saturated and fitted models, multiplied by the scale parameter}
#'      \item{mae}{mean absolute error}
#'      \item{mce}{mean cross entropy, estimated as \eqn{-log-likelihood/N},
#'      where \eqn{N} is the number of observations}
#'      \item{mse}{mean squared error}
#'      \item{r2}{R-squared, calculated as \eqn{1 - deviance/null deviance}}}
#'     }}
#'    \describe{
#'    \item{zero_fit}{a data frame containing the same fit statistics described
#'    for \code{count_fit} but for every possible combination of predictors of
#'    zeroes vs. non-zeroes (\code{n_zero_pred} and \code{zero_pred} replace
#'    \code{n_count_pred} and \code{count_pred}, respectively).}}
#'    \describe{
#'    \item{fit}{a data frame containing the same fit statistics described for
#'    \code{count_fit} but for all combinations of the best model for each
#'    \code{n_count_pred} and the best model for each \code{n_zero_pred} listed
#'    in \code{count_fit} and \code{zero_fit}, respectively.}}
#'    \describe{
#'    \item{cv}{a data frame containing cross-validation statistics for all
#'     combinations of the best models for each \code{n_count_pred} and
#'     \code{n_zero_pred} listed in \code{fit}, except AIC is omitted. Each
#'     metric is followed by its standard error.}}
#'  \describe{
#'    \item{test}{if \code{test_data} is provided, a data frame containing
#'     prediction metrics for the best model for each \code{n_count_pred} and
#'     \code{n_zero_pred} combination listed in \code{fit} as applied to the
#'     \code{test_data}.}}
#'  }}
#' }
#'
#' @export
beset_zeroinfl <- function(form, data, test_data = NULL,
                           family = "poisson", link = "logit",
                           p_count_max = 10, p_zero_max = 10, n_cores = 2,
                           n_folds = 5, n_repeats = 5, seed = 42){
  #==================================================================
  # Check family and link arguments
  #------------------------------------------------------------------
  family <- check_family(family)
  link <- check_link("binomial", link)

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
  if(n_drop > 0)
    warning(paste("Dropping", n_drop, "rows with missing data."),
            immediate. = TRUE)
  # Cache name of the response variable
  response <- names(mf)[1]
  # Extract response vector
  y <- mf[[1]]
  # Insure that a zero-inflated count model is appropriate
  ## Insure that all values are integers
  if(!is.integer(y) || any(y < 0)){
    stop(paste("Response variable does not appear to be a count.",
               "Check that all values are non-negative integers.")
    )
  }
  ## Insure that there are zero values
  if(min(y) != 0){
    stop(paste("Observed lower bound does not equal 0.",
               "Transform your outcome to have a floor of 0,",
               "or, if you do not have a floor effect,",
               "use beset_glm instead.")
    )
  }

  #==================================================================
  # Screen for linear dependencies among predictors
  #------------------------------------------------------------------
  mf <- check_lindep(form, mf)

  #==================================================================
  # Check that number of predictors and cv folds is acceptable
  #------------------------------------------------------------------
  p_count <- min(ncol(mf) - 1, p_count_max)
  n <- nrow(mf)
  alt <- check_cv(n, p_count, FALSE, n_folds)
  if(alt$p < p_count){
    p_count <- alt$p
    warning(paste("'p_count_max' argument too high given sample size",
                  ".\n Reducing maximum subset size for count-model to ",
                  p_count, ".", sep = ""), immediate. = TRUE)
  }
  if(n_folds < alt$folds){
    n_folds <- as.integer(alt$folds)
    warning(paste("'n_folds' argument too low given sample size",
                  "and choice of 'p_count_max'",
                  ".\n  Increasing number of cv folds to ", n_folds, ".",
                  sep = ""), immediate. = TRUE)
  }
  p_zero <- min(ncol(mf) - 1, p_zero_max)
  mf_binom <- mf
  mf_binom[,1] <- mf[,1] > 0
  y_binom <- as.factor(mf_binom[,1])
  n <- min(sum(y_binom == levels(y_binom)[1]),
           sum(y_binom == levels(y_binom)[2]))
  alt <- tryCatch(
    check_cv(n, p_zero, TRUE, n_folds),
    error = function(c){
      c$message <- "Not enough 0 values to model zero-inflation process."
      c$call <- NULL
      stop(c)
    })
  if(alt$p < p_zero){
    p_zero <- alt$p
    warning(paste("'p_zero_max' argument too high given counts of zero values",
                  ".\n Reducing maximum subset size for zero-model to ",
                  p_zero, ".", sep = ""), immediate. = TRUE)
  }
  if(n_folds < alt$folds){
    n_folds <- as.integer(alt$folds)
    warning(paste("'n_folds' argument too low given counts of zero values",
                  "and choice of 'p_zero_max'",
                  ".\n  Increasing number of cv folds to ", n_folds, ".",
                  sep = ""), immediate. = TRUE)
  }
  if(n_folds > alt$folds){
    n_folds <- as.integer(alt$folds)
    warning(paste("'n_folds' argument too high given counts of zero values",
                  ".\n  Decreasing number of cv folds to ", n_folds, ".\n",
                  sep = ""), immediate. = TRUE)
  }
  #======================================================================
  # Make list of all possible formulas with number of predictors <= p_max
  #----------------------------------------------------------------------
  count_pred <- lapply(1:p_count, function(x)
    combn(names(mf)[2:ncol(mf)], x, simplify = FALSE))
  count_pred <- unlist(sapply(count_pred, function(vars){
    sapply(vars, function(x) paste0(x, collapse = " + "))}))
  count_pred <- c("1", count_pred)
  count_form_list <- paste(response, "~", count_pred, "| 1")

  zero_pred <- lapply(1:p_zero, function(x)
    combn(names(mf)[2:ncol(mf)], x, simplify = FALSE))
  zero_pred <- unlist(sapply(zero_pred, function(vars){
    sapply(vars, function(x) paste0(x, collapse = " + "))}))
  zero_pred <- c("1", zero_pred)
  zero_form_list <- paste(response, "~", "1 |", zero_pred)

  #======================================================================
  # Determine number of predictors for each model in list
  #----------------------------------------------------------------------
  n_count_pred <- sapply(count_pred[-1], function(x){
    x <- unlist(strsplit(x, split = " + ", fixed = TRUE))
    sum(sapply(x, function(y){
      if(is.factor(mf[,y]))
        length(levels(mf[,y])) - 1
      else
        1
    }))
  })
  n_count_pred <- c(intercept = 0, n_count_pred)
  n_zero_pred <- sapply(zero_pred[-1], function(x){
    x <- unlist(strsplit(x, split = " + ", fixed = TRUE))
    sum(sapply(x, function(y){
      if(is.factor(mf[,y]))
        length(levels(mf[,y])) - 1
      else
        1
    }))
  })
  n_zero_pred <- c(intercept = 0, n_zero_pred)

  #==================================================================
  # Obtain fit for every model predicting zero vs. non-zero values,
  # with a constant predicting counts.
  #------------------------------------------------------------------
  cl <- parallel::makeCluster(n_cores)
  parallel::clusterExport(cl, c("mf", "family", "link", "test_data",
                                "predict_metrics"), envir=environment())
  zero_fit_stats <- parallel::parLapplyLB(cl, zero_form_list, function(form){
    fit <- try(suppressWarnings(
      pscl::zeroinfl(stats::formula(form), data = mf, dist = family,
                     link = link)),
      silent = TRUE)
    aic <- dev <- mae <- mce <- mse <- r2 <- NA_real_
    if(class(fit) == "zeroinfl"){
      aic <- stats::AIC(fit)
      metrics <- predict_metrics(fit, mf)
      dev <- metrics$deviance
      mae <- metrics$mean_absolute_error
      mce <- metrics$mean_cross_entropy
      mse <- metrics$mean_squared_error
      r2 <- metrics$R_squared
    }
    list(aic = aic, dev = dev, mae = mae, mce = mce, mse = mse, r2 = r2)
  })
  stat_names <- names(zero_fit_stats[[1]])
  get_stat <- function(stat_name, stat_list){
    vapply(stat_list, function(x) x[[stat_name]], 0.0)
  }
  zero_fit_stats <- lapply(stat_names, get_stat, stat_list = zero_fit_stats)
  names(zero_fit_stats) <- stat_names
  zero_fit_stats <- dplyr::bind_cols(
    n_zero_pred = n_zero_pred,
    zero_pred = zero_pred,
    tibble::as_data_frame(zero_fit_stats)
  )
  best_zero_subsets <- dplyr::group_by(zero_fit_stats, n_zero_pred)
  best_zero_subsets <- dplyr::filter(best_zero_subsets,
                                     mce == min(mce, na.rm = TRUE))
  best_zero_subsets <- dplyr::ungroup(best_zero_subsets)

  #=========================================================================
  # Obtain fit for every model predicting counts, with a constant predicting
  # zero vs. non-zero values.
  #-------------------------------------------------------------------------
  count_fit_stats <- parallel::parLapplyLB(cl, count_form_list, function(form){
    fit <- try(suppressWarnings(
      pscl::zeroinfl(stats::formula(form), data = mf, dist = family,
                     link = link)),
      silent = TRUE)
    aic <- dev <- mae <- mce <- mse <- r2 <- NA_real_
    if(class(fit) == "zeroinfl"){
      aic <- stats::AIC(fit)
      metrics <- predict_metrics(fit, mf)
      dev <- metrics$deviance
      mae <- metrics$mean_absolute_error
      mce <- metrics$mean_cross_entropy
      mse <- metrics$mean_squared_error
      r2 <- metrics$R_squared
    }
    list(aic = aic, dev = dev, mae = mae, mce = mce, mse = mse, r2 = r2)
  })
  count_fit_stats <- lapply(stat_names, get_stat, stat_list = count_fit_stats)
  names(count_fit_stats) <- stat_names
  count_fit_stats <- dplyr::bind_cols(
    n_count_pred = n_count_pred,
    count_pred = count_pred,
    tibble::as_data_frame(count_fit_stats)
  )
  best_count_subsets <- dplyr::group_by(count_fit_stats, n_count_pred)
  best_count_subsets <- dplyr::filter(best_count_subsets,
                                     mce == min(mce, na.rm = TRUE))
  best_count_subsets <- dplyr::ungroup(best_count_subsets)

  #======================================================================
  # Make formulas for all combinations of best zero and best count pred
  #----------------------------------------------------------------------
  fit_info <- cbind(
    expand.grid(n_count_pred = best_count_subsets$n_count_pred,
                n_zero_pred = best_zero_subsets$n_zero_pred),
    expand.grid(count_pred = best_count_subsets$count_pred,
                zero_pred = best_zero_subsets$zero_pred,
                stringsAsFactors = FALSE),
    stringsAsFactors = FALSE)
  fit_info$form <- paste(response, "~", fit_info$count_pred,
                             "|", fit_info$zero_pred)
  fit_info <- dplyr::select(fit_info, -count_pred, -zero_pred)
  fit_stats <- parallel::parLapplyLB(cl, fit_info$form, function(form){
    fit <- try(suppressWarnings(
      pscl::zeroinfl(stats::formula(form), data = mf, dist = family,
                     link = link)),
      silent = TRUE)
    aic <- dev <- mae <- mce <- mse <- r2 <- NA_real_
    if(class(fit) == "zeroinfl"){
      aic <- stats::AIC(fit)
      metrics <- predict_metrics(fit, mf)
      dev <- metrics$deviance
      mae <- metrics$mean_absolute_error
      mce <- metrics$mean_cross_entropy
      mse <- metrics$mean_squared_error
      r2 <- metrics$R_squared
    }
    list(aic = aic, dev = dev, mae = mae, mce = mce, mse = mse, r2 = r2)
  })
  fit_stats <- lapply(stat_names, get_stat, stat_list = fit_stats)
  names(fit_stats) <- stat_names
  fit_stats <- dplyr::bind_cols(
    fit_info,
    tibble::as_data_frame(fit_stats)
  )

  #======================================================================
  # Perform cross-validation on best models
  #----------------------------------------------------------------------
  cv_stats <- dplyr::select(fit_stats, n_count_pred:form)
  cv_stats <- dplyr::arrange(cv_stats, n_zero_pred, n_count_pred)
  fold_ids <- create_folds(y, n_folds, n_repeats, seed)
  metrics <- parallel::parLapply(cl, fold_ids, function(i, form_list){
    lapply(form_list, function(form){
      fit <- try(suppressWarnings(
        pscl::zeroinfl(stats::formula(form), data = mf[i,], dist = family,
                       link = link)),
        silent = TRUE)
      if(class(fit) == "zeroinfl"){
        predict_metrics(fit, mf[-i,])
      } else {
        list(deviance = NA_real_,
             mean_absolute_error = NA_real_,
             mean_cross_entropy = NA_real_,
             mean_squared_error = NA_real_,
             R_squared = NA_real_)
      }
    })
  }, form_list = cv_stats$form)
  parallel::stopCluster(cl)
  #======================================================================
  # Derive cross-validation statistics
  #----------------------------------------------------------------------
  metrics <- purrr::transpose(metrics)
  metrics <- purrr::map(metrics, purrr::transpose)
  metrics <- purrr::modify_depth(metrics, 2, purrr::as_vector)
  cv_mean <- purrr::modify_depth(metrics, 2, ~mean.default(.,na.rm = TRUE))
  cv_mean <- purrr::transpose(cv_mean)
  cv_mean <- purrr::map(cv_mean, purrr::as_vector)
  cv_mean <- tibble::as_data_frame(cv_mean)
  names(cv_mean) <- c("dev", "mae", "mce", "mse", "r2")
  cv_se <- purrr::modify_depth(metrics, 2, function(x){
    sd(x, na.rm = TRUE) / sqrt(length(x[!is.na(x)]))
  })
  cv_se <- purrr::transpose(cv_se)
  cv_se <- purrr::map(cv_se, purrr::as_vector)
  cv_se <- tibble::as_data_frame(cv_se)
  names(cv_se) <- c("dev_SE", "mae_SE", "mce_SE", "mse_SE", "r2_SE")
  cv_stats <- dplyr::bind_cols(cv_stats, cv_mean, cv_se)

  #======================================================================
  # Compute prediction statistics for independent test set
  #----------------------------------------------------------------------
  test_stats <- NULL
  if(!is.null(test_data)){
    metrics <- lapply(cv_stats$form, function(form){
      fit <- try(suppressWarnings(
        pscl::zeroinfl(stats::formula(form), data = mf, dist = family,
                       link = link)), silent = TRUE)
      if(class(fit) == "zeroinfl"){
        predict_metrics(fit, test_data)
      } else {
        list(deviance = NA_real_,
             mean_absolute_error = NA_real_,
             mean_cross_entropy = NA_real_,
             mean_squared_error = NA_real_,
             R_squared = NA_real_)
      }
    })
    metrics <- purrr::transpose(metrics)
    metrics <- purrr::map(metrics, purrr::as_vector)
    metrics <- tibble::as_data_frame(metrics)
    names(metrics) <- c("dev", "mae", "mce", "mse", "r2")
    test_stats <- dplyr::bind_cols(
      dplyr::select(cv_stats, n_count_pred, n_zero_pred, form),
                            metrics)
  }
  #======================================================================
  # Construct beset_zeroinfl object
  #----------------------------------------------------------------------
  structure(list(params = list(family = family,
                               link = link,
                               n_folds = n_folds,
                               n_repeats = n_repeats,
                               seed = seed),
                 model_data = mf,
                 stats = list(count_fit = count_fit_stats,
                              zero_fit = zero_fit_stats,
                              fit = fit_stats,
                              cv = cv_stats,
                              test = test_stats)), class = "beset_zeroinfl")
}
