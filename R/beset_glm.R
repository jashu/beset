#' Best Subset Selection for Generalized Linear Models
#'
#' \code{beset_glm} performs best subset selection using repeated
#' cross-validation to find the optimal number of predictors for several
#' families of generalized linear models.
#'
#' \code{beset_glm} performs best subset selection for generalized linear
#' models, fitting a separate model for each possible combination of predictors
#' (all models that contain exactly 1 predictor, all models that contain
#' exactly 2 predictors, and so forth). For each number of predictors,
#' \code{beset_glm} first picks the model with the best fit and then
#' estimates how well this model predicts new data using \code{k}-fold
#' cross-validation (how well, on average, a model trained using \eqn{k - 1}
#' folds predicts the left-out fold).
#'
#' @section Cross-validation details:
#' \code{beset_glm} uses \code{\link{create_folds}} to randomly partition the
#' data set into \code{n_folds} * \code{n_repeats} folds within strata (factor
#' levels for factor outcomes, percentile-based groups for numeric outcomes).
#' This insures that the folds will be matched in terms of the outcome's
#' frequency distribution. \code{beset_glm} also insures the reproducibility of
#' your analysis by requiring a \code{seed} to the random number generator as
#' one of its arguments.
#'
#' @section List of available families and link functions:
#' \describe{
#'  \item{\code{"gaussian"}}{The Gaussian family accepts the links
#'    \code{"identity"} (default), \code{"log"}, and \code{"inverse"}.}
#'  \item{\code{"binomial"}}{The binomial family accepts the links
#'    \code{"logit"} (default), \code{"probit"}, \code{"cauchit"}, \code{"log"}, and
#'    \code{"cloglog"} (complementary log-log).}
#'  \item{\code{"poisson"}}{The Poisson family accepts the links \code{"log"}
#'    (default), \code{"sqrt"}, and \code{"identity"}.}
#'  \item{\code{"negbin"}}{The negative binomial family accepts the links
#'    \code{"log"} (default), \code{"sqrt"}, and \code{"identity"}.}
#'  }
#'
#' @section Warnings:
#' \enumerate{
#'  \item \code{beset_glm} handles missing data by performing listwise deletion.
#'   No other options for handling missing data are provided at this time. The
#'   user is encouraged to deal with missing values prior to running this
#'   function.
#'  \item \code{beset_glm} is intended for use with additive models only.  An
#'  exhaustive search over the space of possible interactions and/or non-linear
#'  effects is computationally prohibitive, but I hope to offer a greedy search
#'  option in the future. In the meantime and in general, I would recommend the
#'    \href{https://cran.r-project.org/web/packages/earth/index.html}{MARS}
#'    technique.
#'  \item \code{beset_glm} is best suited for searching over a small number of
#'  predictors (less than 10). For a large number of predictors (more than 20),
#'  \code{\link{beset_elnet}} is recommended instead. However, note that
#'  \code{\link{beset_elnet}} only works with a more restricted set of
#'  distributions.
#' }
#'
#' @name beset_glm
#' @import parallel
#' @import purrr
#' @importFrom utils combn
#' @import dplyr
#'
#'
#' @seealso \code{\link[stats]{glm}},
#' \code{\link[base]{set.seed}}, \code{\link[MASS]{glm.nb}}
#'
#' @param form A model \code{\link[stats]{formula}}.
#'
#' @param data Either a \code{\link{data_partition}} object containing data sets
#' to be used for both model training and testing, or a single data frame that
#' will be used for model training and cross-validation.
#'
#' @param family Character string naming the error distribution to be used in
#' the model. Available families are listed under 'List of available families
#' and link functions'.
#'
#' @param link (Optional) character string naming the link function to be used in
#' the model. Available links and their defaults differ by \code{family} and are
#' listed under 'List of available families and link functions'.
#'
#' @param p_max Maximum number of predictors to attempt to fit. Default is 10.
#'
#' @param force_in (Optional) vector containing the names or indices of any
#' predictor variables that should be included in every model. (Note that if
#' there is an intercept, it is forced into every model by default.)
#'
#' @param n_folds Integer indicating the number of folds to use for
#' cross-validation.
#'
#' @param n_reps Integer indicating the number of times cross-validation should
#' be repeated (with different randomized fold assignments).
#'
#' @param seed An integer used to seed the random number generator when
#' assigning observations to folds.
#'
#' @param epsilon \code{Numeric} value of positive convergence tolerance ε; the
#' iterations converge when \eqn{|dev - dev_{old}|/(|dev| + 0.1) < ε}. Default
#' is \code{1e-8}.
#'
#' @param maxit \code{Integer} giving the maximal number of IWLS iterations.
#' Default is 25.
#'
#' @param skinny \code{Logical} value indicating whether or not to return a
#' "skinny" model. If \code{FALSE} (the default), the return object will include
#' a copy of the model \code{\link[stats]{terms}}, \code{data},
#' \code{contrasts}, and a record of the \code{xlevels} of the factors used in
#' fitting. If these features are not needed, setting \code{skinny = TRUE} will
#' prevent these copies from being made.
#'
#' @param n_cores Integer value indicating the number of workers to run in
#' parallel during subset search and cross-validation. By default, this will
#' be set to one fewer than the maximum number of physical cores you have
#' available, as indicated by \code{\link[parallel]{detectCores}}. Set to 1 to
#' disable parallel processing.
#'
#' @param parallel_type (Optional) character string indicating the type of
#' parallel operation to be used, either \code{"fork"} or \code{"sock"}. If
#' omitted and \code{n_cores > 1}, the default is \code{"sock"} for Windows and
#' \code{"fork"} for any other OS.
#'
#' @param cl (Optional) \code{\link[parallel]{parallel}} or
#' \code{\link[snow]{snow}} cluster for use if \code{parallel_type = "sock"}.
#' If not supplied, a cluster on the local machine is automatically created.
#'
#' @inheritParams stats::glm
#'
#' @return A "beset_glm" object with the following components:
#' \describe{
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
#'      \item{test}{if \code{test_data} is provided, a data frame
#'      containing prediction metrics for the best model for each \code{n_pred}
#'      listed in \code{fit} as applied to the \code{test_data}.}
#'      }
#'    }
#'   \item{fold_assignments}{list giving the row indices for the holdout
#'    observations for each fold and/or repetition of cross-validation}
#'   \item{n_folds}{number of folds used in cross-validation}
#'   \item{n_reps}{number of repetitions used in cross-validation}
#'   \item{family}{name of error distribution used in the model}
#'   \item{link}{name of link function used in the model}
#'   \item{terms}{the \code{\link[stats]{terms}} object used}
#'   \item{data}{the \code{data} argument}
#'   \item{offset}{the offset vector used}
#'   \item{contrasts}{(where relevant) the contrasts used}
#'   \item{xlevels}{(where relevant) a record of the levels of the factors used
#'        in fitting}
#'     }
NULL

#' @rdname beset_glm
#' @export
beset_glm <- function(form, data, family = "gaussian", link = NULL,
                      p_max = 10, force_in = NULL,
                      nest_cv = FALSE, n_folds = 10, n_reps = 10, seed = 42,
                      contrasts = NULL, offset = NULL, weights = NULL,
                      start = NULL, etastart = NULL, mustart = NULL,
                      epsilon = 1e-8, maxit = 25, skinny = FALSE,
                      n_cores = NULL, parallel_type = NULL, cl = NULL){

  #==================================================================
  # Check family argument and identify appropriate model fit function
  #------------------------------------------------------------------
  family <- check_family(family)
  fitter <- if(family == "negbin") "glm_nb" else "glm.fit"

  #==================================================================
  # Check for missing values and linear dependence
  #------------------------------------------------------------------
  if(inherits(data, "data.frame")){
    data <- model.frame(form, data = data)
    n_omit <- length(attr(data, "na.action"))
    if(n_omit > 0){
      warning(paste("Dropping", n_omit, "rows with missing data."),
              immediate. = TRUE)
      attr(data, "na.action") <- NULL
    }
    data <- check_lindep(data)
    mf <- model.frame(form, data = data)
    terms <- terms(mf)
    xlevels = .getXlevels(terms, mf)

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
  } else if(is.null(n_cores) || n_cores > 1){
      if(nest_cv && is.null(parallel_type)) parallel_type <- "sock"
      parallel_control <- setup_parallel(
        parallel_type = parallel_type, n_cores = n_cores, cl = cl)
      n_cores <- parallel_control$n_cores
      cl <- parallel_control$cl
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
          beset_glm(form = form, data = x, family = family, link = link,
                    p_max = p_max, force_in = force_in, contrasts = contrasts,
                    nest_cv = FALSE, n_folds = n_folds, n_reps = 1, seed = seed,
                    start = start, etastart = etastart, mustart = mustart,
                    epsilon = epsilon, maxit = maxit, n_cores = 1)
        })
      } else {
        parallel::parLapply(cl, all_partitions, function(x){
          beset_glm(form = form, data = x, family = family, link = link,
                    p_max = p_max, force_in = force_in, contrasts = contrasts,
                    nest_cv = FALSE, n_folds = n_folds, n_reps = 1, seed = seed,
                    start = start, etastart = etastart, mustart = mustart,
                    epsilon = epsilon, maxit = maxit, n_cores = 1)
        })
      }
    } else {
      lapply(all_partitions, function(x){
        beset_glm(form = form, data = x, family = family, link = link,
                  p_max = p_max, force_in = force_in, contrasts = contrasts,
                  nest_cv = FALSE, n_folds = n_folds, n_reps = 1, seed = seed,
                  start = start, etastart = etastart, mustart = mustart,
                  epsilon = epsilon, maxit = maxit, n_cores = 1)
      })
    }
    if(skinny){
      terms <- data <- contrasts <- xlevels <- NULL
    }
    out <- structure(
      list(
        beset = nested_cv, fold_assignments = fold_ids,
        n_folds = n_folds, n_reps = n_reps,
        family = family, link = link, terms = terms, data = data,
        offset = offset, contrasts = contrasts, xlevels = xlevels
      ),
      class = c("nested", "beset", "glm")
    )
    if(!is.null(cl)) stopCluster(cl)
    return(out)
  }

  #==================================================================
  # Create list of arguments for model
  #------------------------------------------------------------------
  m <- beset:::make_args(form = form, data = data, family = family, link = link,
                 contrasts = contrasts, weights = weights, offset = offset,
                 start = start, etastart = etastart, mustart = mustart,
                 epsilon = epsilon, maxit = maxit)

  #==================================================================
  # Set number of cross-validation folds and reps
  #------------------------------------------------------------------
  n_obs <- length(m$train$y)
  cv_params <- set_cv_par(n_obs = n_obs, n_folds = n_folds, n_reps = n_reps)
  n_reps <- cv_params$n_reps; n_folds <- cv_params$n_folds

  #======================================================================
  # Get all subsets (see utils_subset.R)
  #----------------------------------------------------------------------
  all_subsets <- get_subsets(m, force_in, p_max)

  #======================================================================
  # Obtain fit statistics for every model
  #----------------------------------------------------------------------
  train_test_stats <- if (n_cores > 1L) {
    if(is.null(cl)){
      parallel::mclapply(all_subsets$pred_idx, get_subset_stats, m = m,
                         fitter= fitter, mc.cores = n_cores)
    } else {
      parallel::parLapply(cl, all_subsets$pred_idx, get_subset_stats, m = m,
                          fitter = fitter)
    }
  } else lapply(all_subsets$pred_idx, get_subset_stats, m = m, fitter = fitter)
  fit_stats <- transpose(train_test_stats)$fit_stats %>%
    transpose %>% simplify_all %>% as_data_frame %>%
    mutate_all(function(x) case_when(abs(x) < 1e-6 ~ 0, TRUE ~ x))
  fit_stats <- dplyr::bind_cols(all_subsets, dplyr::as_data_frame(fit_stats))
  test_stats <- NULL
  test_preds <- NULL
  if(!is.null(train_test_stats[[1]]$test_stats)){
    test_stats <- transpose(train_test_stats)$test_stats %>%
      transpose %>% simplify_all %>% as_data_frame %>%
      mutate_all(function(x) case_when(abs(x) < 1e-6 ~ 0, TRUE ~ x))
    test_stats <- bind_cols(all_subsets, as_data_frame(test_stats))
    test_preds <- transpose(train_test_stats)$test_preds
  }
  m$test$yhat <- test_preds

  #======================================================================
  # Obtain model with best fit for each number of parameters
  #----------------------------------------------------------------------
  cv_stats <- group_by(fit_stats, n_pred) %>%
    filter(rsq == max(rsq)) %>%
    ungroup() %>%
    select(pred_idx:n_pred) %>%
    arrange(n_pred)
  best_fits <- lapply(cv_stats$pred_idx, function(j){
    model <- c(list(x = m$train$x[, j, drop = FALSE]), m$train[-1])
    out <- do.call(fitter, model)
    out <- c(out, model[c("x", "offset", "control", "intercept")])
    class(out) <- c("glm", "lm")
    if(family == "negbin") class(out) <- c("negbin", class(out))
    out
  })

  #======================================================================
  # Perform cross-validation on best models
  #----------------------------------------------------------------------
  cv_results <- if (n_cores > 1L) {
    if(is.null(cl)){
      parallel::mclapply(best_fits, validate, n_folds = n_folds,
                         n_reps = n_reps, seed = seed, mc.cores = n_cores)
    } else {
      parallel::parLapply(cl, best_fits, validate, n_folds = n_folds,
                          n_reps = n_reps, seed = seed)
    }
  } else map(best_fits, ~ validate(., n_folds = n_folds, n_reps = n_reps,
                                       seed = seed)
  )
  cv_results <- cv_results %>% map("stats") %>% transpose %>% as_data_frame

  cv_stats <- bind_cols(cv_stats, cv_results)


  #======================================================================
  # Construct beset_glm object
  #----------------------------------------------------------------------
  stats <- list(fit = fit_stats, cv = cv_stats, test = test_stats)
  parameters <- m$train
  if(skinny){
    terms <- data <- contrasts <- xlevels <- fold_ids <- NULL
  } else {
    fold_ids <- create_folds(m$train$y, n_folds, n_reps, seed)
  }
  if(!is.null(cl)) stopCluster(cl)
  structure(
    list(stats = stats,
         parameters = parameters,
         fold_assignments = fold_ids, n_folds = n_folds, n_reps = n_reps,
         family = family, link = link, terms = terms, data = data,
         offset = offset, contrasts = contrasts, xlevels = xlevels),
    class = c("beset", "glm")
  )
}

#' @export
#' @rdname beset_glm
beset_lm <- function(form, data, p_max = 10, force_in = NULL,
                     weights = NULL, contrasts = NULL, offset = NULL,
                     nest_cv = FALSE, n_folds = 10, n_reps = 10,  seed = 42,
                     n_cores = NULL, parallel_type = NULL, cl = NULL){
  do.call(beset_glm, as.list(match.call())[-1])
}

