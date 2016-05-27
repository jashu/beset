#' Beset GLM with Elasticnet Regularization
#'
#' \code{beset_elnet} is a wrapper to \code{\link[glmnet]{glmnet}} for fitting
#'  generalized linear modesl via penalized maximum likelihood, providing
#'  automated data preprocessing and selection of both the elasticnet penalty
#'  and regularization parameter through cross-validation.
#'
#' @param form A formula of the form y ~ x1 + x2 + ...
#'
#' @param train_data Data frame containing the training data set and all of the
#' variables specified in \code{form}.
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
#' @param seed Seed for random number generator.
#'
#' @param n_cores Number of cores to use for parallel execution. If not
#' specified, the number of cores is set to the value of
#' \code{options("cores")}, if specified, or to approximately half the number of
#' cores detected by the \code{parallel} package. To determine the theoretical
#' maximum number of cores you have available, see
#' \code{\link[parallel]{detectCores}}, but note that the actual number of cores
#' available may be less. See \href{
#' https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf}{
#' \code{parallel} package vignette} for more information.
#'
#' @seealso \code{\link[glmnet]{glmnet}}, \code{\link[caret]{train}}
#'
#' @export

beset_elnet <- function(form, train_data, standardize = TRUE,
                        method = c("best", "1SE"),
                        max_features = 1000, alpha = c(.05, .5, .95),
                        n_folds = 5, n_repeats = 5, seed = 42, n_cores = NULL){
  #===================================================================
  # Register parallel backend
  #-------------------------------------------------------------------
  if(is.null(n_cores)) n_cores <- parallel::detectCores() %/% 2
  doParallel::registerDoParallel(cores = n_cores)
  #==================================================================
  # ETL x and y from formula and data frame
  #------------------------------------------------------------------
  if("data.frame" %in% class(train_data)){
    mf <- model.frame(form, data = train_data)
    x <- as.matrix(mf[,2:ncol(mf)])
    y <- mf[,1]
  } else {
    x <- train_data[, 2:ncol(train_data)]
    y <- train_data[, 1]
  }
  if(any(is.na(y))){
    warning("Cases missing response will be deleted.", immediate. = TRUE)
    x <- x[!is.na(y),]
    y <- na.omit(y)
  }
  if(standardize){
    x <- apply(x, 2, scale)
  }

  #==================================================================
  # For very high-dimensional data, implement recursively
  #------------------------------------------------------------------
  while(ncol(x) > max_features){
    n_features <- ncol(x)
    vars_to_keep <- NULL
    for(start in seq(1, n_features, max_features)){
      ptm <- proc.time()
      end <- min(start + max_features - 1, n_features)
      temp <- cbind(y, x[, start:end])
      colnames(temp)[1] <- all.vars(form)[1]
      temp_model <- beset_elnet(form, train_data = temp,
                                max_features = max_features,
                                standardize = FALSE, method = "best")
      temp_summary <- summary(temp_model)
      vars_to_keep <- c(vars_to_keep, temp_summary$var_imp$variable[
        abs(temp_summary$var_imp$best) > 0])
      elapsed_time <- (proc.time() - ptm)[3]
      runs_remain <- (n_features - end) / max_features
      time_remain <- elapsed_time * runs_remain
      hrs_remain <- time_remain %/% (60 * 60)
      time_remain <- time_remain %% (60 * 60)
      min_remain <- time_remain %/% 60
      message(paste(end, "vars screened;",
                    length(vars_to_keep), "vars selected."))
      message(paste("Estimated time remaining:",
                    hrs_remain, "hours,",
                    min_remain, "minutes."))
      save(vars_to_keep, file = ".temp.RData")
    }
    x <- x[, vars_to_keep]
    if(ncol(x) > max_features) max_features <- 2 * max_features
  }

  #==================================================================
  # Assign cross-validation folds
  #------------------------------------------------------------------
  folds <- matrix(nrow = length(y), ncol = n_repeats)
  set.seed(seed, kind = "default")
  for(r in 1:n_repeats){
    folds[, r] <- caret::createFolds(y, k = n_folds, list = FALSE)
  }
  #======================================================================
  # Train glmnet
  #----------------------------------------------------------------------
  family <- "gaussian"
  if(dplyr::n_distinct(y) == 2){
    family <- "binomial"
  }
  results <- data.frame()
  pb <- txtProgressBar(min = 0, max = length(alpha) * n_repeats, style = 3)
  i <- 1
  for(a in seq_along(alpha)){
    for(r in 1:n_repeats){
      temp <- glmnet::cv.glmnet(x, y, foldid = folds[, r], parallel = TRUE,
                        alpha = alpha[a], family = family, standardize = FALSE)
      temp <- cbind(alpha = alpha[a], as.data.frame(temp[1:2]))
      results <- rbind(results, temp)
      setTxtProgressBar(pb, i)
      i <- i + 1
    }
  }

  #======================================================================
  # Extract results and determine largest alpha and smallest lambda within 1 SE
  #----------------------------------------------------------------------
  results <- dplyr::mutate(results, lambda = round(lambda, 3))
  results <- dplyr::group_by(results, alpha, lambda)
  results <- dplyr::summarise(results, cve = mean(cvm),
                              cve_se = sqrt(var(cvm) / length(cvm)))
  results <- na.omit(results)
  results <- dplyr::ungroup(results)
  results <- dplyr::mutate(results, cve_lo = cve - cve_se,
                           cve_hi = cve + cve_se)
  best_result <- dplyr::filter(results, cve == min(cve))
  best_results <- dplyr::filter(results, cve < best_result$cve_hi)
  best_result_1SE <- dplyr::filter(best_results, alpha == max(alpha))
  best_result_1SE <- dplyr::filter(best_result_1SE, lambda == max(lambda))
  best_alpha <- best_result$alpha
  best_lambda <- best_result$lambda
  best_alpha_1SE <- best_result_1SE$alpha
  best_lambda_1SE <- best_result_1SE$lambda

  #======================================================================
  # Fit best model and best model within 1 SD
  #----------------------------------------------------------------------
  best_model <- NULL; best_model_1SE <- NULL
  if("best" %in% method){
    best_model <- glmnet::glmnet(x, y, family = family,
                                 alpha = best_alpha,
                                 standardize = FALSE)
  }
  if("1SE" %in% method){
    if(best_alpha == best_alpha_1SE){
      best_model_1SE <- best_model
      } else {
        best_model_1SE <- glmnet::glmnet(x, y, family = family,
                                         alpha = best_alpha_1SE,
                                         standardize = FALSE)
      }
  }
  structure(list(results = results,
                 best_alpha = best_alpha,
                 best_alpha_1SE = best_alpha_1SE,
                 best_lambda = best_lambda,
                 best_lambda_1SE = best_lambda_1SE,
                 best_model = best_model,
                 best_model_1SE = best_model_1SE),
            class = "beset_elnet")
}
