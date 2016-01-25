#' Beset GLM with Elasticnet Regularization
#'
#' \code{beset_glmnet} is a wrapper to \code{\link[glmnet]{glmnet}} via
#' \code{\link[caret]{train}} for fitting generalized linear modesl via
#' penalized maximum likelihood, providing automated data preprocessing and
#' selection of both the elasticnet penalty and regularization parameter through
#' cross-validation.
#'
#' \code{beset_glmnet} provides automatic centering and scaling of the predictor
#' variables (which is a must for regularized regression) and optional multiple
#' imputation of missing values. Imputation is repeated within each
#' cross-validation fold and thus variance among imputations becomes a part of
#' the cross-validation error.
#'
#' @param form A formula of the form y ~ x1 + x2 + ...
#'
#' @param train_data Data frame containing the training data set and all of the
#' variables specified in \code{form}.
#'
#' @param test_data Optional data frame containing an independent test data set
#' and all of the variables specified in \code{form}.
#'
#' @param alpha Numeric vector of alpha values between 0 and 1 to use as tuning
#' parameters. \code{alpha = 0} results in ridge regression, and \code{alpha =
#' 1} results in lasso regression. Values in between result in a mixture of L1
#' and L2 penalties. (Values closer to 0 weight the L2 penalty more heavily,
#' and values closer to 1 weight the L1 penalty more heavily.)
#'
#' @param cv_control Integer vector of length two indicating the number of
#' cross-validation folds followed by the number of times cross-validation
#' should be repeated. Default is \code{c(5,5)}, i.e., 5-fold cross-validation,
#' repeated 5 times.
#'
#' @param impute Optional character string giving a method for imputing missing
#' values. This must be one of the strings "none" (the default, which will
#' result in listwise deletion of any cases with missing values), "knn"
#' (K-nearest neighbors), "bag" (bagged trees), or "median".
#'
#' @param seed Seed for random number generator.
#'
#' @param cores Number of cores to use for parallel execution. If not specified,
#' computations will be performed serially using a single core. To determine
#' the theoretical maximum number of cores you have available, see
#' \code{\link[parallel]{detectCores}}, but note that the actual number of cores
#' available may be less. See \href{
#' https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf}{
#' \code{parallel} package vignette} for more information.
#'
#' @seealso \code{\link[glmnet]{glmnet}}, \code{\link[caret]{train}}
#'
#' @export

beset_glmnet <- function(form, train_data, test_data = NULL,
                         alpha = c(.05, .5, .95),
                         cv_control = c(5,5),
                         impute = "none",
                         seed = 42, parallel = FALSE){
  #===================================================================
  # Register parallel processing
  #-------------------------------------------------------------------
  if(parallel){
    doMC::registerDoMC()
  }
  #===================================================================
  # Set up preprocessing
  #-------------------------------------------------------------------
  preprocess <- c("center", "scale")
  impute_options <- c("none", "knn", "bag", "median")
  if(!impute %in% impute_options){
    stop(paste("Unrecognized method of imputation. Impute must be one of",
               "\"none\", \"knn\", \"bag\", or \"median\"."))
  }
  if(impute != "none"){
    preprocess <- c(preprocess, paste(impute, "Impute", sep = ""))
  }
  #==================================================================
  # ETL x and y from formula and data frame
  #------------------------------------------------------------------
  mf <- model.frame(form, data = train_data,
                    na.action = ifelse(impute == "none", "na.omit", "na.pass"))
  x <- as.matrix(mf[,2:ncol(mf)])
  y <- mf[,1]
  if(any(is.na(y))){
    warning("Cases missing response will be deleted.", immediate. = TRUE)
    x <- x[!is.na(y),]
    y <- na.omit(y)
  }
  if(length(unique(y)) == 2 && !is.factor(y)){
    y <- binary_to_factor(y)
  }
  #==================================================================
  # Make summary function
  #------------------------------------------------------------------
  glmnet_summary <- function(data, lev = NULL, model = NULL){
    if(is.factor(y)){
      return(c(caret::defaultSummary(data, lev, model),
               caret::twoClassSummary(data, lev, model),
               caret::mnLogLoss(data, lev, model)))
    } else {
      return(caret::defaultSummary(data, lev, model))
    }
  }
  #======================================================================
  # Initialize tuning grid
  #----------------------------------------------------------------------
  glmnet_grid <- expand.grid(alpha = alpha,
                             lambda = c(seq(.001, .009, .001),
                                        seq(.01, .09, .01),
                                        seq(.1, 1, .1)))
  # =====================================================================
  # Make random seeds to insure reproducibility when cross-validation
  # is performed in parallel across independent cores.
  #----------------------------------------------------------------------
  set.seed(seed)
  seeds <- vector(mode = "list", length = cv_control[1] * cv_control[2] + 1)
  for(i in 1:(length(seeds) - 1)){
    seeds[[i]] <- sample.int(length(seeds)^2, length(alpha))
  }
  seeds[[length(seeds)]] <- seed
  #======================================================================
  # Set up train control
  #----------------------------------------------------------------------
  ctrl <- caret::trainControl(method = "repeatedcv",
                                      number = cv_control[1],
                                      repeats = cv_control[2],
                                      verboseIter = TRUE,
                                      classProbs = is.factor(y),
                                      summaryFunction = glmnet_summary,
                                      seeds = seeds)
  #======================================================================
  # Train glmnet
  #----------------------------------------------------------------------
  metric = ifelse(is.factor(y), "logLoss", "RMSE")
  glmnet_train <- caret::train(x = x,
                               y = y,
                               method = "glmnet",
                               preProcess = preprocess,
                               metric = metric,
                               trControl = ctrl,
                               tuneGrid = glmnet_grid)
  #======================================================================
  # Extract results, transform standard deviations to standard errors,
  # and determine largest alpha and smallest lambda within 1 SE.
  #----------------------------------------------------------------------
  results <- glmnet_train$results
  SD_cols <- grep("SD$", names(results))
  results[, SD_cols] <- results[, SD_cols] / sqrt(cv_control[1] * cv_control[2])
  names(results)[SD_cols] <- gsub("SD$", "SE", names(results)[SD_cols])
  best_result <- results[which.min(results[, metric]),
                         c(metric, paste(metric, "SE", sep = ""))]
  oneSE_results <- results[results[, metric] < sum(best_result),]
  oneSE_results <-
    oneSE_results[oneSE_results$alpha == max(oneSE_results$alpha),]
  oneSE_results <-
    oneSE_results[oneSE_results$lambda == min(oneSE_results$lambda),]
  #======================================================================
  # Fit 1SE model.
  #----------------------------------------------------------------------
  best_model <- glmnet_train$finalModel
  best_model_1SE <- best_model
  if(oneSE_results$alpha != best_model$tuneValue$alpha){
    glmnet_grid <- data.frame(alpha = oneSE_results$alpha,
                              lambda = oneSE_results$lambda)
    ctrl$method <- "none"
    glmnet_train <- caret::train(x = x,
                                 y = y,
                                 method = "glmnet",
                                 preProcess = preprocess,
                                 metric = metric,
                                 trControl = ctrl,
                                 tuneGrid = glmnet_grid)
    best_model_1SE <- glmnet_train$finalModel
  }
  best_model_1SE$lambdaOpt <- oneSE_results$lambda
  best_model_1SE$tuneValue$lambda <- oneSE_results$lambda
  structure(list(results = results,
                 best_model = best_model,
                 best_model_1SE = best_model_1SE),
            class = "beset_glmnet")
}
