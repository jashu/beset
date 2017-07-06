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
#' \code{beset_glm} uses \code{\link[caret]{createMultiFolds}} to randomly
#' partition the data set into \code{n_folds} * \code{n_repeats} folds within
#' strata (factor levels for factor outcomes, percentile-based groups for
#' numeric outcomes). This insures that the folds will be matched in terms of
#' the outcome's frequency distribution. \code{beset_glm} also insures the
#' reproducibility of your analysis by requiring a \code{seed} to the random
#' number generator as one of its arguments.
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
#'  \item \code{beset_glm} is intended for use with additive models only. If you
#'   include interaction terms in the model formula passed to \code{form}, they
#'   will be ignored (e.g. \code{A * B * C} will be treated as \code{A + B + C}
#'   and \code{A + B + A:B} will be treated as \code{A + B}). An exhaustive
#'   search over the space of possible interactions and/or non-linear effects is
#'   computationally prohibitive, but I hope to offer a greedy search option in
#'   the future. In the meantime and in general, I would recommend the
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
#' @importFrom purrr at_depth
#' @importFrom purrr map
#' @importFrom purrr map_chr
#' @importFrom purrr map_int
#' @importFrom purrr flatten_dbl
#' @importFrom purrr reduce
#' @importFrom purrr transpose
#' @importFrom utils combn
#' @import dplyr
#'
#'
#' @seealso \code{\link[caret]{createFolds}}, \code{\link[stats]{glm}},
#' \code{\link[base]{set.seed}}, \code{\link[MASS]{glm.nb}}
#'
#' @param form A model \code{\link[stats]{formula}}.
#'
#' @param data Data frame with the variables in \code{form} and the data
#' to be used for model fitting.
#'
#' @param test_data Optional data frame with the variables in \code{form} and
#' the data to be used for model validation.
#'
#' @param family Character string naming the error distribution to be used in
#' the model. Available families are listed under 'List of available families
#' and link functions'.
#'
#' @param link Optional character string naming the link function to be used in
#' the model. Available links and their defaults differ by \code{family} and are
#' listed under 'List of available families and link functions'.
#'
#' @param p_max Maximum number of predictors to attempt to fit.
#'
#' @param n_cores Integer value indicating the number of workers to run in
#' parallel during subset search and cross-validation. By default, this will
#' be set to 2. To determine the theoretical maximum number of cores you have
#' available, see \code{\link[parallel]{detectCores}}, but note that the actual
#' number of cores available may be less. See
#' \code{\link[parallel]{parallel-package}} for more information.
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
#' @return A "beset_glm" object with the following components:
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
NULL

#' @rdname beset_glm
#' @export
beset_glm <- function(form, data, test_data = NULL, p_max = 10,
                      family = "gaussian", link = NULL,
                      n_cores = 2, n_folds = 10, n_repeats = 10, seed = 42){

  #==================================================================
  # Check family argument and set up link function if specified
  #------------------------------------------------------------------
  family <- check_family(family)
  link <- if(!is.null(link)){
    check_link(family, link)
  } else {
    switch(family,
           binomial = "logit",
           gaussian = "identity",
           poisson = "log",
           negbin = "log")
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
    warning(paste("Dropping", n_drop, "rows with missing data."),
            immediate. = TRUE)
  }
  # cache name of the response variable
  response <- names(mf)[1]
  # extract response vector
  y <- mf[[1]]
  # insure y is a factor if family is binomial
  if(family == "binomial" && !is.factor(y)) y <- factor(y)

  #==================================================================
  # Screen for linear dependencies among predictors
  #------------------------------------------------------------------
  mf <- check_lindep(form, mf)

  #==================================================================
  # Check that number of predictors and cv folds is acceptable
  #------------------------------------------------------------------
  p <- min(ncol(mf) - 1L, p_max)
  n <- if(family == "binomial"){
    min(sum(y == levels(y)[1]), sum(y == levels(y)[2]))
  } else {
    nrow(mf)
  }
  alt <- check_cv(n, p, family == "binomial", n_folds)
  if(alt$p < p){
    p <- alt$p
    warning(paste("'p_max' argument too high for your sample size",
                  ".\n  Reducing maximum subset size to ", p, ".",
                  sep = ""), immediate. = TRUE)
  }
  if(n_folds < alt$folds){
    n_folds <- as.integer(alt$folds)
    warning(paste("'n_folds' argument too low for your sample size ",
                  "and choice of 'p_max'",
                  ".\n  Increasing number of cv folds to ", n_folds, ".",
                  sep = ""), immediate. = TRUE)
  }
  #======================================================================
  # Make list of all possible formulas with number of predictors <= p_max
  #----------------------------------------------------------------------
  pred <- lapply(1:p, function(x)
    combn(names(mf)[2:ncol(mf)], x, simplify = FALSE))
  pred <- reduce(pred, c)
  form_list <- map_chr(pred, ~ paste0(.x, collapse = " + "))
  form_list <- c("1", form_list)
  form_list <- paste(response, "~", form_list)

  #======================================================================
  # Determine number of predictors for each model in list
  #----------------------------------------------------------------------
  get_npred <- function(pred_name){
    pred_data <- mf[[pred_name]]
    if(is.factor(pred_data)) length(levels(pred_data)) - 1L else 1L
  }
  n_pred <- pred %>% at_depth(2, get_npred) %>% map_int(reduce, `+`)
  n_pred <- c(0L, n_pred)

  #======================================================================
  # Eliminate formulae that exceed p_max
  #----------------------------------------------------------------------
  form_list <- form_list[n_pred <= p]
  n_pred <- n_pred[n_pred <= p]

  #======================================================================
  # Obtain fit statistics for every model
  #----------------------------------------------------------------------
  cl <- parallel::makeCluster(n_cores)
  parallel::clusterExport(cl, c("mf", "family", "link", "test_data", "glm_nb",
                                "fit_glm"), envir = environment())
  fit_stats <- parallel::parLapplyLB(cl, form_list, function(form){
    fit <- try(fit_glm(mf, form, family, link), silent = TRUE)
    aic <- dev <- mae <- mce <- mse <- r2 <- NA_real_
    if(!inherits(fit, "try-error")){
       aic <- stats::AIC(fit)
       dev <- stats::deviance(fit)
       mae <- mean(abs(stats::residuals(fit, type = "response")))
       mce <- -1 * as.numeric(stats::logLik(fit)) / nrow(mf)
       mse <- mean(stats::residuals(fit, type = "response")^2)
       r2 <- 1 - fit$deviance / fit$null.deviance
    }
    list(aic = aic, dev = dev, mae = mae, mce = mce, mse = mse, r2 = r2)
    }) %>% transpose() %>% map(flatten_dbl) %>% as_data_frame() %>%
    mutate_all(function(x) case_when(abs(x) < 1e-6 ~ 0, TRUE ~ x))
  stat_names <- names(fit_stats)
  fit_stats <- bind_cols(
    n_pred = n_pred,
    form = form_list,
    as_data_frame(fit_stats)
    )

  #======================================================================
  # Obtain model with best fit for each number of parameters
  #----------------------------------------------------------------------
  cv_stats <- group_by(fit_stats, n_pred) %>%
    filter(r2 == max(r2)) %>%
    ungroup() %>%
    select(n_pred, form) %>%
    arrange(n_pred)

  #======================================================================
  # Perform cross-validation on best models
  #----------------------------------------------------------------------
  set.seed(seed)
  fold_ids <- caret::createMultiFolds(y, k = n_folds, times = n_repeats)
  metrics <- parallel::parLapplyLB(cl, fold_ids, function(i, form_list){
    lapply(form_list, function(form){
      fit <- fit_glm(mf[i,], form, family, link)
      stats <- try(predict_metrics(fit, test_data = mf[-i,]), silent = FALSE)
      if(inherits(stats, "prediction_metrics")){
        stats
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

  cv_mean <-  map(1:nrow(cv_stats), function(i)
    metrics %>% map(i) %>% transpose() %>% map(flatten_dbl) %>%
      map(mean, na.rm = TRUE)
  ) %>% transpose() %>% map(flatten_dbl) %>% as_data_frame()
  names(cv_mean) <- stat_names[-1]

  se <- function(x){
    sd(x, na.rm = TRUE) / sqrt(length(x[!is.na(x)]))
  }

  cv_se <-  map(1:nrow(cv_stats), function(i)
    metrics %>% map(i) %>% transpose() %>% map(flatten_dbl) %>% map(se)
  ) %>% transpose() %>% map(flatten_dbl) %>% as_data_frame() %>%
    mutate_all(function(x) case_when(abs(x) < 1e-6 ~ 0, TRUE ~ x))
  names(cv_se) <- paste(stat_names[-1], "SE", sep = "_")

  cv_stats <- bind_cols(cv_stats, cv_mean, cv_se)

  #======================================================================
  # Compute prediction statistics for independent test set
  #----------------------------------------------------------------------
  test_stats <- NULL
  if(!is.null(test_data)){
    metrics <- lapply(cv_stats$form, function(form){
      fit <- fit_glm(test_data, form, family, link)
      stats <- tryCatch(
        predict_metrics(fit, test_data = test_data),
        error = function(c){
          c$message <- paste("unable to compute 'test_data' prediction metrics",
                             "for one or more subsets because\n", c$message)
          c$call <- NULL
          warning(c)
          })
      if(class(stats) == "prediction_metrics") stats
      else list(deviance = NA_real_,
                mean_absolute_error = NA_real_,
                mean_cross_entropy = NA_real_,
                mean_squared_error = NA_real_,
                R_squared = NA_real_)
      })
    metrics <- transpose(metrics) %>%
      at_depth(1, as_vector) %>%
      as.data.frame()
    names(metrics) <- c("dev", "mae", "mce", "mse", "r2")
    test_stats <- bind_cols(select(cv_stats, n_pred, form), metrics)
  }
  #======================================================================
  # Construct beset_glm object
  #----------------------------------------------------------------------
  structure(list(params = list(family = family,
                               link = link,
                               n_folds = n_folds,
                               n_repeats = n_repeats,
                               seed = seed),
                 model_data = mf,
                 stats = list(fit = fit_stats,
                              cv = cv_stats,
                              test = test_stats)),
            class = "beset_glm")
}

#' @export
#' @rdname beset_glm
beset_lm <- function(form, data, test_data = NULL, p_max = 10,
                     n_cores = 2, n_folds = 10, n_repeats = 10,  seed = 42){
  beset_glm(form, data, test_data = test_data, p_max = p_max,
            n_cores = n_cores, n_folds = n_folds, n_repeats = n_repeats,
            seed = seed)
}

