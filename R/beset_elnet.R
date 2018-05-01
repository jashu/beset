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
#' @param alpha_seq \code{Numeric} vector of alpha values between 0 and 1 to
#' use as tuning parameters. \code{alpha = 0} results in ridge regression, and
#' \code{alpha = 1} results in lasso regression. Values in between result in a
#' mixture of L1 and L2 penalties. (Values closer to 0 weight the L2 penalty
#' more heavily, and values closer to 1 weight the L1 penalty more heavily.) The
#' default is to try three alpha values: 0.01 (emphasis toward ridge penalty),
#' 0.99 (emphasis toward lasso penalty), and 0.5 (equal mixture of L1 and L2).
#'
#' @param lambda_seq (Optional) \code{numeric} vector giving a decreasing
#' sequence of \eqn{lambda} values. If omitted, the program computes its own
#' \code{lambda_seq} based on \code{n_lambda} and \code{lambda_min_ratio}.
#'
#' @param n_lambda Number of lambdas to be used in a search. Defaults to
#' \code{100}.
#'
#' @param lambda_min_ratio (Optional) minimum lambda used in lambda search,
#' specified as a ratio of lambda_max (the smallest lambda that drives all
#' coefficients to zero). Default if omitted: if the number of observations is
#' greater than the number of variables, then \code{lambda_min_ratio} is set to
#' 0.0001; if the number of observations is less than the number of variables,
#' then \code{lambda_min_ratio} is set to 0.01.
#'
#' @param offset A vector of length nobs that is included in the linear
#' predictor. Useful for the "poisson" family (e.g. log of exposure time), or
#' for refining a model by starting at a current fit. Default is NULL.
#'
#' @param weights (Optional) \code{numeric} vector of observation weights
#' of length equal to the number of cases.
#'
#' @param epsilon Convergence threshold for coordinate descent. Each inner
#' coordinate-descent loop continues until the maximum change in the objective
#' after any coefficient update is less than thresh times the null deviance.
#' Defaults value is \code{1e-7}.
#'
#' @param maxit Maximum number of passes over the data for all lambda values;
#' default is \code{10^5}.
#'
#' @param standardize Logical flag for x variable standardization, prior to
#' fitting the model sequence. The coefficients are always returned on the
#' original scale. Default is \code{standardize = TRUE}. If variables are in the
#' same units already, you might not wish to standardize.
#'
#' @param contrasts {Optional} \code{list}. See the \code{contrasts.arg} of
#' \code{\link[stats]{model.matrix.default}}.
#'
#' @return A "beset_elnet" or "nested_elnet" object with the following
#' components:
#' \enumerate{
#'  \item\describe{
#'    \item{params}{list of parameters used in function call that will be
#'    needed to reproduce results}
#'    }
#'   \item\describe{
#'     \item{model_data}{data frame containing training data used to identify
#'      best subsets.}
#'  }
#'  \item\describe{
#'    \item{stats}{a list with three data frames:
#'      \describe{
#'        \item{fit}{statistics for every possible combination of predictors:
#'          \describe{
#'            \item{n_pred}{the total number of predictors in model; note that
#'               the number of predictors for a factor variable corresponds to the
#'               number of factor levels minus 1}
#'            \item{form}{formula for model}
#'            \item{aic}{\eqn{-2*log-likelihood + k*npar}, where \eqn{npar}
#'              represents the number of parameters in the fitted model, and
#'              \eqn{k = 2}}
#'            \item{dev}{twice the difference between the log-likelihoods of the
#'              saturated and fitted models, multiplied by the scale parameter}
#'            \item{mae}{mean absolute error}
#'            \item{mce}{mean cross entropy, estimated as
#'              \eqn{-log-likelihood/N}, where \eqn{N} is the number of
#'              observations}
#'            \item{mse}{mean squared error}
#'            \item{r2}{R-squared, calculated as
#'              \eqn{1 - deviance/null deviance}}
#'            }
#'          }
#'      \item{cv}{a data frame containing cross-validation statistics
#'      for the best model for each \code{n_pred} listed in \code{fit_stats}.
#'      Each metric is computed using \code{\link{predict_metrics}}, with
#'      models fit to \eqn{n-1} folds and predictions made on the left-out fold.
#'      Each metric is followed by its standard error. The data frame
#'      is otherwise the same as that documented for \code{fit}, except
#'      AIC is omitted.}
#'      \item{test_stats}{if \code{test_data} is provided, a data frame
#'      containing prediction metrics for the best model for each \code{n_pred}
#'      listed in \code{fit} as applied to the \code{test_data}.}
#'       }
#'     }
#'   }
#' }
#'
#' @inheritParams beset_glm
#' @seealso \code{\link[glmnet]{glmnet}}
#'
#' @import purrr
#' @export

beset_elnet <- function(form, data, family = "gaussian",
                        alpha_seq = c(.01, .5, .99), lambda_seq = NULL,
                        nest_cv = FALSE, n_folds = 10, n_reps = 10, seed = 42,
                        weights = NULL, offset = NULL, contrasts = NULL,
                        n_lambda = 100, lambda_min_ratio = NULL, p_max = NULL,
                        standardize = TRUE, epsilon = 1e-7, maxit = 10^5,
                        parallel_type = NULL, n_cores = NULL, cl = NULL){

  #==================================================================
  # Check family argument
  #------------------------------------------------------------------
  family <- tryCatch(
    match.arg(family, c("binomial", "gaussian", "poisson")),
    error = function(c){
      c$message <- gsub("arg", "family", c$message)
      c$call <- NULL
      stop(c)
    })

  #======================================================================
  # Set up parallel operations
  #----------------------------------------------------------------------
  if(is.null(n_cores) || n_cores > 1){
    parallel_control <- setup_parallel(
      parallel_type = parallel_type, n_cores = n_cores, cl = cl)
    have_mc <- parallel_control$have_mc
    n_cores <- parallel_control$n_cores
    cl <- parallel_control$cl
  }

  #======================================================================
  # Recursive function for performing nested cross-validation
  #----------------------------------------------------------------------
  if(nest_cv){
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
    fold_ids <- create_folds(data[[y]], n_folds, n_reps, seed)
    all_partitions <- lapply(fold_ids, function(i){
      suppressWarnings(
        data_partition(train = data[-i,], test = data[i,], y = y, x = x,
                       offset = o, weights = w)
      )
    })
    out <- if(n_cores > 1L){
      if(have_mc){
        parallel::mclapply(all_partitions, function(x){
          beset_elnet(form = form, data = x, family = family,
                      contrasts = contrasts, alpha_seq = alpha_seq,
                      lambda_seq = lambda_seq, n_lambda = n_lambda,
                      lambda_min_ratio = lambda_min_ratio,
                      standardize = standardize, epsilon = epsilon,
                      p_max = p_max, maxit = maxit, n_folds = n_folds,
                      n_reps = 1, seed = seed, n_cores = 1)
          })
        } else {
          parallel::parLapply(cl, all_partitions, function(x){
            beset_elnet(form = form, data = x, family = family,
                        contrasts = contrasts, alpha_seq = alpha_seq,
                        lambda_seq = lambda_seq, n_lambda = n_lambda,
                        lambda_min_ratio = lambda_min_ratio,
                        standardize = standardize, epsilon = epsilon,
                        p_max = p_max, maxit = maxit, n_folds = n_folds,
                        n_reps = 1, seed = seed, n_cores = 1)
          })
        }
      } else {
        lapply(all_partitions, function(x){
        beset_elnet(form = form, data = x, family = family,
                    contrasts = contrasts, alpha_seq = alpha_seq,
                    lambda_seq = lambda_seq, n_lambda = n_lambda,
                    lambda_min_ratio = lambda_min_ratio,
                    standardize = standardize, epsilon = epsilon,
                    p_max = p_max, maxit = maxit, n_folds = n_folds,
                    n_reps = 1, seed = seed, n_cores = 1)
      })
      }
    attr(out, "fold_assignments") <- fold_ids
    class(out) <- "nested_elnet"
    return(out)
  }


  #==================================================================
  # Create model matrices
  #------------------------------------------------------------------
  m <- map(alpha_seq, ~ make_args(
    form = form, data = data, family = family, link = NULL, nlambda = n_lambda,
    contrasts = contrasts, weights = weights, offset = offset,
    epsilon = epsilon, maxit = maxit, alpha = .x, lambda = lambda_seq,
    standardize = standardize)
  )
  n_obs <- nrow(m[[1]]$train$x)
  names(m) <- alpha_seq
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
    if(!is.null(p_max)){
      m[[i]]$train$dfmax <<- p_max
      if(length(m[[i]]$test)) m[[i]]$test$dfmax <<- p_max
    }
  })

  #======================================================================
  # Obtain fit statistics
  #----------------------------------------------------------------------
  fits <- if (n_cores > 1L) {
    if (have_mc) {
      parallel::mclapply(m, function(x) do.call(glmnet::glmnet, x$train),
                         mc.cores = min(n_cores, length(m)))
    } else {
      parallel::parLapply(cl, m, function(x) do.call(glmnet::glmnet, x$train))
    }
  } else lapply(m, function(x) do.call(glmnet::glmnet, x$train))
  names(fits) <- alpha_seq
  lambda_seq <- map(fits, "lambda")
  alpha_seq <- map2(alpha_seq, lambda_seq, ~ rep(.x, length(.y)))
  fit_stats <- data.frame(alpha = unlist(alpha_seq),
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
    test_stats <- data.frame(alpha = unlist(alpha_seq),
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
    test_preds <- y_hats
  }

  #==================================================================
  # Obtain cross-validation stats
  #------------------------------------------------------------------
  cv_params <- set_cv_par(n_obs, n_folds, n_reps)
  n_folds <- cv_params$n_folds; n_reps <- cv_params$n_reps
  cv_results <- map(fits, ~ validate(.x, n_folds = n_folds, n_reps = n_reps,
                                     seed = seed)) %>%
    map("stats") %>% transpose %>% map(reduce, c) %>% as_data_frame
  cv_stats <- fit_stats %>% select(alpha, lambda) %>% bind_cols(cv_results)

  #======================================================================
  # Construct beset_elnet object
  #----------------------------------------------------------------------
  stats <- list(fit = fit_stats, cv = cv_stats, test = test_stats)
  fit_parameters <- m[[1]]$train
  fit_parameters$alpha <- NULL
  test_parameters <- m[[1]]$test
  test_parameters$predictions <- test_preds
  cv_parameters <- list(n_folds = n_folds, n_reps = n_reps, seed = seed,
                        n_cores = n_cores)
  structure(list(stats = stats,
                 parameters = list(fit = fit_parameters,
                                   cv = cv_parameters,
                                   test = test_parameters)
                 ), class = "beset_elnet")
}
