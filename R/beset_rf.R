#' Beset Random Forest
#'
#' \code{beset_rf} is a wrapper to \code{\link[randomForest]{randomForest}} that
#' estimates predictive performance of the random forest using repeated k-fold
#' cross-validation. \code{beset_rf} insures that the correct arguments are
#' provided to \code{\link[randomForest]{randomForest}} and that enough
#' information is retained for compatibility with \code{beset} methods such as
#' variable \code{\link{importance}} and partial \code{\link{dependence}}.
#'
#' @param data Data frame containing the variables in the model.
#'
#' @param n_trees Number of trees. Defaults to 500.
#'
#' @param sample_rate Row sample rate per tree (from \code{0 to 1}). Defaults to
#' \code{0.6320000291}.
#'
#' @param mtry (Optional) \code{integer} number of variables randomly sampled
#' as candidates at each split. If omitted, defaults to the square root of the
#' number of predictors for classification and one-third the number of
#' predictors for regression.
#'
#' @param min_obs_in_node (Optional) \code{integer} number specifying the
#' fewest allowed observations in a terminal node. If omitted, defaults to 1 for
#' classification and 5 for regression.
#'
#' @param class_wt Priors of the classes. Ignored for regression.
#'
#' @param x A \code{"beset_rf"} object to plot
#'
#' @param metric Prediction metric to plot. Options are mean squared error
#' (\code{"mse"}) or R-squared (\code{"rsq"}) for regression, and
#' misclassification error (\code{"err.rate"}) for classification. Default
#' \code{"auto"} plots MSE for regression and error rate for classification.
#'
#' @inheritParams beset_glm
#' @inheritParams randomForest::randomForest
#'
#' @return A "beset_rf" object with the following components:
#' \describe{
#'   \item{forests}{list of "randomForest" objects for each fold and repetition}
#'   \item{stats}{a "cross_valid" object giving cross-validation metrics}
#'   \item{data}{the data frame used to train random forest}
#'  }
#'
#' @examples
#' data("prostate", package = "beset")
#' rf <- beset_rf(tumor ~ ., data = prostate)
#' @export
beset_rf <- function(form, data, n_trees = 500, sample_rate = 0.6320000291,
                     mtry = NULL, min_obs_in_node = NULL,
                     n_folds = 10, n_reps = 10, seed = 42,
                     class_wt = NULL, cutoff = NULL, strata = NULL,
                     parallel_type = NULL, n_cores = NULL, cl = NULL){
  #======================================================================
  # Set up x and y arguments for randomForest.default
  #----------------------------------------------------------------------
  data <- mutate_if(data, is.logical, factor)
  data <- model.frame(form, data = data)
  n_omit <- length(attr(data, "na.action"))
  if(n_omit > 0){
    warning(paste("Dropping", n_omit, "rows with missing data."),
            immediate. = TRUE)
    attr(data, "na.action") <- NULL
  }
  x <- data[-1]
  y <- data[[1]]
  names(y) <- row.names(data)
  if(is.factor(y)){
    if(n_distinct(y) > 2){
      stop(
        paste("Multinomial classification not supported by beset at this time.",
              "Please use the `randomForest` package directly for this.",
              sep = "\n")
      )
    }
  } else {
    if(n_distinct(y) == 2) y <- factor(y)
  }
  #======================================================================
  # Set up parallel operations
  #----------------------------------------------------------------------
  if(!is.null(cl)){
    if(!inherits(cl, "cluster")) stop("Not a valid parallel socket cluster")
    n_cores <- length(cl)
  } else if(is.null(n_cores) || n_cores > 1){
    if(is.null(parallel_type)) parallel_type <- "sock"
    parallel_control <- setup_parallel(
      parallel_type = parallel_type, n_cores = n_cores, cl = cl)
    n_cores <- parallel_control$n_cores
    cl <- parallel_control$cl
  }
  #======================================================================
  # Set up cross-validation
  #----------------------------------------------------------------------
  n_obs <- length(y)
  cv_par <- beset:::set_cv_par(n_obs, n_folds, n_reps)
  n_folds <- cv_par$n_folds; n_reps <- cv_par$n_reps
  fold_ids <- beset:::create_folds(y = y, n_folds = n_folds, n_reps = n_reps,
                                   seed = seed)
  #======================================================================
  # Set up arguments for randomForest
  #----------------------------------------------------------------------
  rf_par <- list(
    ntree = n_trees,
    mtry = if(!is.null(y) && !is.factor(y)){
      max(floor(ncol(x)/3), 1)
    } else {
      floor(sqrt(ncol(x)))
    },
    replace = FALSE, classwt = class_wt,
    cutoff = if(is.null(cutoff) && is.factor(y)){
      n_class <- length(levels(y))
      rep(1/n_class, n_class)
    } else cutoff, strata = strata,
    sampsize = ceiling(sample_rate * nrow(x)),
    nodesize = if(is.null(min_obs_in_node)){
      if(!is.null(y) && !is.factor(y)) 5 else 1
    } else min_obs_in_node,
    importance = TRUE, localImp = TRUE, keep.forest = TRUE
  )
  train_data <- lapply(fold_ids, function(i)
    c(list(x = x[-i, , drop = FALSE],
           y = y[-i],
           xtest = x[i, , drop = FALSE],
           ytest = y[i]),
      rf_par)
  )
  #======================================================================
  # Train and cross-validate random forests
  #----------------------------------------------------------------------
  cv_fits <- if(n_cores > 1L){
    if(is.null(cl)){
      parallel::mclapply(train_data, function(x){
        set.seed(seed); a <- do.call(randomForest::randomForest, x)
      }, mc.cores = n_cores)
    } else {
      parallel::parLapply(cl, train_data, function(x){
        set.seed(seed); do.call(randomForest::randomForest, x)
      })
    }
  } else {
    lapply(train_data, function(x){
      set.seed(seed); do.call(randomForest::randomForest, x)
    })
  }
  if(!is.null(cl)) parallel::stopCluster(cl)
  type <- if(is.factor(y)) "prob" else "response"
  y_hat <- map2(cv_fits, train_data,
                ~ as.matrix(predict(.x, .y$xtest, type = type)))
  if(is.factor(y)){
    y_hat <- map(y_hat, ~ .x[, 2, drop = FALSE])
    family <- "binomial"
  } else {
    family <- "gaussian"
  }
  cv_stats <- beset:::get_cv_stats(y = y, y_hat = y_hat, family = family,
                           n_folds = n_folds, n_reps = n_reps)
  fold_assignments <- get_fold_ids(fold_ids, n_reps)
  cv_stats <- structure(
    c(cv_stats, list(
      fold_assignments = fold_assignments,
      parameters = list(family = family,
                        n_obs = n_obs,
                        n_folds = n_folds,
                        n_reps = n_reps,
                        seed = seed,
                        y = y))),
    class = "cross_valid"
  )
  structure(
    list(
      forests = cv_fits,
      stats = cv_stats,
      data = data
    ), class = c("beset_rf", "beset")
  )
}

#' @export
#' @describeIn beset_rf Plot OOB and holdout MSE, R-squared, or error rate as a
#'  function of number of trees in forest
#'
plot.beset_rf <- function(x, metric = c("auto", "mse", "rsq", "err.rate"), ...){
  metric <- tryCatch(
    match.arg(metric, c("auto", "mse", "rsq", "err.rate")),
    error = function(c){
      c$message <- gsub("arg", "metric", c$message)
      c$call <- NULL
      stop(c)
    }
  )
  if(metric == "auto"){
    metric <- if(x$forests[[1]]$type == "regression") "mse" else "err.rate"
  }
  if(metric %in% c("mse", "rsq") && x$forests[[1]]$type != "regression"){
    warning(paste(metric, " plots not available for classification.\n",
                  "Misclassification rate plotted instead.", sep = ""))
    metric <- "err.rate"
  }
  oob <- map(x$forests, ~ .x[[metric]])
  if(inherits(oob[[1]], "matrix")) oob <- map(oob, ~ .x[,1])
  oob <- oob %>% transpose %>% simplify_all %>% map_dbl(mean)
  oob <- data_frame(
    sample = "Out-of-Bag", n_trees = seq(1:length(oob)), mean = oob
  )
  cv <- map(x$forests, ~ .x$test[[metric]])
  if(inherits(cv[[1]], "matrix")) cv <- map(cv, ~ .x[,1])
  cv <- cv %>% transpose %>% simplify_all %>% map_dbl(mean)
  cv <- data_frame(
    sample = "Test Holdout", n_trees = seq(1:length(cv)), mean = cv
  )
  data <- bind_rows(oob, cv)
  y_lab <- switch(metric,
                  err.rate = ylab("Misclassification Rate"),
                  mse = ylab("Mean Squared Error"),
                  rsq = ylab(bquote(~R^2)))
  p <- ggplot(data = data, aes(x = n_trees, y = mean, color = sample)) +
    theme_classic() + xlab("Number of trees") + y_lab +
    geom_line(size = 1) + theme(legend.title = element_blank())
  suppressWarnings(print(p))
}
