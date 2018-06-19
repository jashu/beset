#' Beset GLM with Elasticnet Regularization
#'
#' \code{beset_elnet} is a wrapper to \code{\link[glmnet]{glmnet}} for fitting
#'  generalized linear models via penalized maximum likelihood, providing
#'  automated data preprocessing and selection of both the elastic-net penalty
#'  and regularization parameter through repeated k-fold cross-validation.
#'
#' @param form A model \code{\link[stats]{formula}}.
#'
#' @param data Either a \code{\link{data_partition}} object containing data sets
#' to be used for both model training and testing, or a single data frame that
#' will be used for model training only.
#'
#' @param family \code{Character} string naming the error distribution to be
#' used in the model. Currently supported options are \code{"gaussian"}
#' (default), \code{"binomial"}, and \code{"poisson"}.
#'
#' @param skinny \code{Logical} value indicating whether or not to return a
#' "skinny" model. If \code{FALSE} (the default), the return object will include
#' a copy of the model \code{\link[stats]{terms}}, \code{data},
#' \code{contrasts}, and a record of the \code{xlevels} of the factors used in
#' fitting. This information will be necessary if you apply
#' \code{\link{predict.beset_elnet}} to new data. If this feature is not needed,
#' setting \code{skinny = TRUE} will prevent these copies from being made.
#'
#' @param nest_cv \code{Logical} value indicating whether or not to perform a
#' nested cross-validation that isolates the cross-validation used for tuning
#' \code{alpha} and \code{lambda} from the cross-validation used to estimate
#' prediction error. Setting to \code{TRUE} will increase run time considerably
#' (by a factor equal to the number of folds), but useful for estimating
#' uncertainty in the tuning procedure. Defaults to \code{FALSE}.
#'
#' @param alpha \code{Numeric} vector of alpha values between 0 and 1 to
#' use as tuning parameters. \code{alpha = 0} results in ridge regression, and
#' \code{alpha = 1} results in lasso regression. Values in between result in a
#' mixture of L1 and L2 penalties. (Values closer to 0 weight the L2 penalty
#' more heavily, and values closer to 1 weight the L1 penalty more heavily.) The
#' default is to try three alpha values: 0.01 (emphasis toward ridge penalty),
#' 0.99 (emphasis toward lasso penalty), and 0.5 (equal mixture of L1 and L2).
#'
#' @param n_lambda Number of lambdas to be used in a search. Defaults to
#' \code{100}.
#'
#' @param lambda_min_ratio (Optional) minimum \eqn{lambda} used in \eqn{lambda}
#' search, specified as a ratio of \code{lambda_max} (the smallest \eqn{lambda}
#' that drives all coefficients to zero). Default if omitted: if the number of
#' observations is greater than the number of variables, then
#' \code{lambda_min_ratio} is set to 0.0001; if the number of observations is
#' less than the number of variables, then \code{lambda_min_ratio} is set to
#' 0.01.
#'
#' @param offset (Optional) vector of length equal to the number of observations
#' that is included in the linear predictor. Useful for the "poisson" family
#' (e.g. log of exposure time), or for refining a model by starting at a current
#' fit.
#'
#' @param weights (Optional) \code{numeric} vector of observation weights
#' of length equal to the number of cases.
#'
#' @param epsilon Convergence threshold for coordinate descent.
#'
#' @param maxit Maximum number of passes over the data for all lambda values
#'
#' @param standardize Logical flag for x variable standardization, prior to
#' fitting the model sequence. The coefficients are always returned on the
#' original scale. Default is \code{standardize = TRUE}. If variables are in the
#' same units already, you might not wish to standardize.
#'
#' @param contrasts {Optional} \code{list}. See the \code{contrasts.arg} of
#' \code{\link[stats]{model.matrix.default}}.
#'
#' @param remove_collinear_columns \code{Logical}. In case of linearly dependent
#' columns, remove some of the dependent columns. Defaults to FALSE.
#'
#' @return A "beset_elnet" or "nested" object inheriting class "beset_elnet"
#' with the following components:
#'
#' \describe{
#' \item{For "beset_elnet" objects:}{
#' \describe{
#'    \item{stats}{a list with three data frames:
#'      \describe{
#'        \item{fit}{
#'          \describe{
#'            \item{alpha}{value of L1-L2 mixing parameter}
#'            \item{lambda}{value of shrinkage parameter}
#'            \item{auc}{area under curve (binomial models only)}
#'            \item{mae}{mean absolute error (not given for binomial models)}
#'            \item{mce}{mean cross entropy, estimated as
#'              \eqn{-log-likelihood/N}, where \eqn{N} is the number of
#'              observations}
#'            \item{mse}{mean squared error}
#'            \item{rsq}{R-squared, calculated as
#'              \eqn{1 - deviance/null deviance}}
#'            }
#'          }
#'      \item{cv}{a data frame containing cross-validation statistics for each
#'        \code{alpha} and \code{lambda} listed in \code{fit}. If run with
#'        \code{nest_cv = TRUE}, this will correspond to the inner
#'        cross-validation used to select \code{alpha} and \code{lambda}. Each
#'        metric consists of the following list:
#'        \describe{
#'          \item{mean}{mean of the metric calculated on the aggregate holdout
#'            folds for each repetition and averaged across repetitions}
#'          \item{btwn_fold_se}{the variability between all holdout folds, given
#'            as a standard error}
#'          \item{btwn_rep_range}{after aggregating over all hold-out folds
#'            within each repetition, the variability between repetitions, given
#'            as a min-max range}
#'        }
#'      }
#'      \item{test}{if a \code{\link{data_partition}} is provided, or if run
#'      with \code{nest_cv = TRUE}, a data frame containing prediction metrics
#'      for each \code{alpha} and \code{lambda} listed in \code{fit} as applied
#'      to the independent test data or outer cross-validation holdout data}
#'       }
#'     }
#'   \item{glmnet_parameters}{a list of all parameters that were passed to
#'     \code{\link[glmnet]{glmnet}}}
#'  }}}
#'
#'  \describe{
#'  \item{For "nested" objects:}{
#'  \describe{
#'    \item{beset_elnet}{a list of "beset_elnet" objects, one for each train-
#'      test partition of the outer cross-validation procedure, each consisting
#'      of all of the elements listed above}
#'  }}}
#'
#'  \describe{
#'  \item{For both "nested" and unnested "beset_elnet" objects:}{
#'  \describe{
#'    \item{fold_assignments}{list giving the row indices for the holdout
#'    observations for each fold and/or repetition of cross-validation}
#'    \item{n_folds}{number of folds used in cross-validation}
#'    \item{n_reps}{number of repetitions used in cross-validation}
#'    \item{family}{names of error distribution used in the model}
#'    \item{terms}{the \code{\link[stats]{terms}} object used}
#'    \item{data}{the \code{data} argument}
#'    \item{offset}{the offset vector used}
#'    \item{contrasts}{(where relevant) the contrasts used}
#'    \item{xlevels}{(where relevant) a record of the levels of the factors used
#'      in fitting}
#'  }}}
#'
#' @inheritParams beset_glm
#' @seealso \code{\link[glmnet]{glmnet}}
#'
#' @import purrr
#' @export

beset_elnet <- function(
  form, data, family = "gaussian", alpha = c(.01, .5, .99), n_lambda = 100,
  nest_cv = FALSE, n_folds = 10, n_reps = 10, seed = 42,
  remove_collinear_columns = FALSE, skinny = FALSE, standardize = TRUE,
  epsilon = 1e-7, maxit = 10^5, lambda_min_ratio = NULL,
  contrasts = NULL, offset = NULL, weights = NULL,
  parallel_type = NULL, n_cores = NULL, cl = NULL
  ){

  #==================================================================
  # Check variable names and family
  #------------------------------------------------------------------
  check_names(names(data))
  family <- tryCatch(
    match.arg(family, c("binomial", "gaussian", "poisson")),
    error = function(c){
      c$message <- gsub("arg", "family", c$message)
      c$call <- NULL
      stop(c)
    })
  #==================================================================
  # Check for missing values and linear dependence
  #------------------------------------------------------------------
  if(inherits(data, "data.frame")){
    mf <- model.frame(form, data = data)
    n_omit <- nrow(data) - nrow(mf)
    if(n_omit > 0){
      warning(paste("Dropping", n_omit, "rows with missing data."),
              immediate. = TRUE)
      attr(data, "na.action") <- NULL
    }
    terms <- terms(mf)
    xlevels = .getXlevels(terms, mf)
    attr(mf, "terms") <- NULL
    data <- mf
    if(remove_collinear_columns){
      data <- check_lindep(mf)
      mf <- model.frame(form, data = data)
      terms <- terms(mf)
      xlevels = .getXlevels(terms, mf)
    }
  } else if(inherits(data, "data_partition")){
    terms <- terms(data$train)
    xlevels <- .getXlevels(terms, data$train)
  } else {
    stop("`data` argument must inherit class 'data.frame' or 'data_partition'")
  }

  #======================================================================
  # Set up parallel operations
  #----------------------------------------------------------------------
  if(!is.null(cl)){
    if(!inherits(cl, "cluster")) stop("Not a valid parallel socket cluster")
    n_cores <- length(cl)
  } else {
    if(is.null(n_cores) || n_cores > 1){
      if(!nest_cv){
        if(is.null(n_cores)){
          n_cores <- length(alpha)
        } else {
          n_cores <- min(n_cores, length(alpha))
        }
      } else {
        if(is.null(parallel_type)) parallel_type <- "sock"
      }
      parallel_control <- setup_parallel(
        parallel_type = parallel_type, n_cores = n_cores, cl = cl)
      n_cores <- parallel_control$n_cores
      cl <- parallel_control$cl
    }
  }

  #======================================================================
  # Recursive function for performing nested cross-validation
  #----------------------------------------------------------------------
  if(nest_cv){
    if(inherits(data, "data_partition")){
      tryCatch(
        stop(paste("Ambiguous call: do you want to use the test partition found",
                   "in `data`` to estimate test performance, or do you want to",
                   "use nested cross-validation? Either set `nest_cv` to FALSE",
                   "or pass a single data frame instead of a train/test split.",
                   sep = "\n")),
        error = function(c){
          c$call <- NULL
          stop(c)
        })
    }
    o <- NULL; w <- NULL
    if(!is.null(offset)){
      data$offset = offset
      o <- "offset"
    }
    if(!is.null(weights)){
      data$weights = weights
      w <- "weights"
    }
    y <- all.vars(form)[1]
    x <- all.vars(form)[-1]
    if(length(x) == 1 && x == ".") x <- NULL
    n_obs <- nrow(na.omit(data))
    cv_params <- set_cv_par(n_obs, n_folds, n_reps)
    n_folds <- cv_params$n_folds; n_reps <- cv_params$n_reps
    fold_ids <- create_folds(data[[y]], n_folds, n_reps, seed)
    all_partitions <- lapply(fold_ids, function(i){
      suppressWarnings(
        data_partition(train = data[-i,], test = data[i,], y = y, x = x,
                       offset = o, weights = w)
      )
    })
    nested_cv <- if(n_cores > 1L){
      if(is.null(cl)){
        parallel::mclapply(all_partitions, function(x){
          beset_elnet(form = form, data = x, family = family,
                      contrasts = contrasts, alpha = alpha,
                      n_lambda = n_lambda,
                      lambda_min_ratio = lambda_min_ratio, skinny = TRUE,
                      standardize = standardize, epsilon = epsilon,
                      maxit = maxit, n_folds = n_folds, n_reps = 1,
                      seed = seed, n_cores = 1)
          })
        } else {
          parallel::parLapply(cl, all_partitions, function(x){
            beset_elnet(form = form, data = x, family = family,
                        contrasts = contrasts, alpha = alpha,
                        n_lambda = n_lambda,
                        lambda_min_ratio = lambda_min_ratio, skinny = TRUE,
                        standardize = standardize, epsilon = epsilon,
                        maxit = maxit, n_folds = n_folds, n_reps = 1,
                        seed = seed, n_cores = 1)
          })
        }
      } else {
        lapply(all_partitions, function(x){
        beset_elnet(form = form, data = x, family = family,
                    contrasts = contrasts, alpha = alpha,
                    n_lambda = n_lambda,
                    lambda_min_ratio = lambda_min_ratio, skinny = TRUE,
                    standardize = standardize, epsilon = epsilon,
                    maxit = maxit, n_folds = n_folds, n_reps = 1,
                    seed = seed, n_cores = 1)
      })
      }
    if(skinny){
      terms <- data <- contrasts <- xlevels <- NULL
    }
    out <- structure(
      list(
        beset = nested_cv, fold_assignments = fold_ids,
        n_folds = n_folds, n_reps = n_reps,
        family = family, terms = terms, data = data, offset = offset,
        contrasts = contrasts, xlevels = xlevels
      ),
      class = c("nested", "beset", "elnet")
    )
    if(!is.null(cl)) stopCluster(cl)
    return(out)
  }

  #==================================================================
  # Create model matrices
  #------------------------------------------------------------------
  m <- map(alpha, ~ make_args(
    form = form, data = data, family = family, link = NULL, nlambda = n_lambda,
    contrasts = contrasts, weights = weights, offset = offset,
    epsilon = epsilon, maxit = maxit, alpha = .x, standardize = standardize)
  )
  n_obs <- nrow(m[[1]]$train$x)
  cv_params <- set_cv_par(n_obs, n_folds, n_reps)
  names(m) <- alpha
  # glmnet handles intercept differently from other predictors
  # so should not be included in model matrix
  walk(seq_along(m), function(i){
      if(m[[i]]$train$intercept){
        m[[i]]$train$x <<- m[[i]]$train$x[,-1]
        if(length(m[[i]]$test)) m[[i]]$test$x <<- m[[i]]$test$x[,-1]
    }
    if(!is.null(lambda_min_ratio)){
      m[[i]]$train$lambda.min.ratio <<- lambda_min_ratio
      if(length(m[[i]]$test)) m[[i]]$test$lambda.min.ratio <<- lambda_min_ratio
    }
  })

  #======================================================================
  # Obtain fit statistics
  #----------------------------------------------------------------------
  fits <- if (n_cores > 1L) {
    if(is.null(cl)){
      parallel::mclapply(m, function(x) do.call(glmnet::glmnet, x$train),
                         mc.cores = n_cores)
    } else {
      parallel::parLapply(cl, m, function(x) do.call(glmnet::glmnet, x$train))
    }
  } else lapply(m, function(x) do.call(glmnet::glmnet, x$train))
  names(fits) <- alpha
  lambda_seq <- map(fits, "lambda")
  alpha <- map2(alpha, lambda_seq, ~ rep(.x, length(.y)))
  fit_stats <- data.frame(alpha = unlist(alpha),
                          lambda = unlist(lambda_seq))
  y_hats <- map2(fits, m,
                 ~ predict(.x, .y$train$x, type = "response",
                           newoffset = .y$train$offset) %>%
                   as_data_frame) %>% reduce(c)
  y_obs <- if(is.factor(m[[1]]$train$y))
    as.integer(m[[1]]$train$y) - 1 else m[[1]]$train$y
  fit_stats <- bind_cols(
    fit_stats, map(y_hats, ~ predict_metrics_(y_obs, ., family)) %>%
      transpose %>% simplify_all %>% as_data_frame)

  #==================================================================
  # Obtain independent test stats
  #------------------------------------------------------------------
  test_stats <- NULL; test_preds <- NULL
  if(!is.null(m[[1]]$test)){
    test_stats <- data.frame(alpha = unlist(alpha),
                             lambda = unlist(lambda_seq))
    y_hats <- map2(fits, m,
                   ~ predict(.x, .y$test$x, type = "response",
                             newoffset = .y$test$offset) %>%
                     as_data_frame) %>% reduce(c)
    y_obs <- if(is.factor(m[[1]]$test$y))
      as.integer(m[[1]]$test$y) - 1 else m[[1]]$test$y
    test_stats <- bind_cols(
      test_stats, map(y_hats, ~ predict_metrics_(y_obs, ., family)) %>%
        transpose %>% simplify_all %>% as_data_frame)
  }

  #==================================================================
  # Obtain cross-validation stats
  #------------------------------------------------------------------
  cv_results <- map(
    fits, ~ validate(.x, n_folds = n_folds, n_reps = n_reps, seed = seed)) %>%
    map("stats") %>% transpose %>% map(reduce, c) %>% as_data_frame
  cv_stats <- fit_stats %>% select(alpha, lambda) %>% bind_cols(cv_results)

  #======================================================================
  # Construct beset_elnet object
  #----------------------------------------------------------------------
  stats <- list(fit = fit_stats, cv = cv_stats, test = test_stats)
  parameters <- m[[1]]$train
  parameters$alpha <- NULL
  if(skinny){
    terms <- data <- contrasts <- xlevels <- fold_ids <- NULL
  } else {
    fold_ids <- create_folds(y_obs, n_folds, n_reps, seed)
  }
  if(!is.null(cl)) stopCluster(cl)
  structure(
    list(stats = stats,
         parameters = parameters,
         fold_assignments = fold_ids, n_folds = n_folds, n_reps = n_reps,
         family = family, terms = terms, data = data, offset = offset,
         contrasts = contrasts, xlevels = xlevels),
    class = c("beset", "elnet")
  )
}
