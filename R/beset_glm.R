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
#' \code{beset_glm} picks the model with the maximum likelihood and then
#' estimates how well this model predicts new data using \code{k}-fold
#' cross-validation (how well a model trained using \eqn{k - 1} folds
#' predicts the left-out fold).
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
#'  \item \code{beset_glm} is intended for use with additive models only.
#'    While there is no prohibition against the inclusion of interaction or
#'    polynomial terms, this practice is strongly discouraged. At best, this
#'    will result in an inefficient search because \code{beset_glm} performs an
#'    exhaustive search over all possible variable subsets, including subsets
#'    that are hierarchically incomplete, i.e., subsets that contain an
#'    interaction term but are missing one or more of the subterms that comprise
#'    it. At worst, it may return one of these hierarchically incomplete models
#'    as the best model, an undesirable result if one cares about
#'    interpretability. If one wishes the model search to include interaction
#'    and/or non-linear effects, the
#'    \href{https://cran.r-project.org/web/packages/earth/index.html}{MARS}
#'    technique is recommended instead.
#'  \item \code{beset_glm} is best suited for searching over a small number of
#'  predictors (less than 10). For a large number of predictors (more than 20),
#'  \code{\link{beset_elnet}} is recommended instead. However, note that
#'  \code{\link{beset_elnet}} only works with a more restricted set of
#'  distributions.
#' }
#'
#' @name beset_glm
#' @importFrom utils combn
#' @importFrom purrr at_depth
#' @importFrom purrr transpose
#' @importFrom purrr as_vector
#' @import dplyr
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
  family <- tryCatch(match.arg(family, c("binomial", "gaussian", "poisson",
                                         "negbin")),
                     error = function(c){
                       c$message <- gsub("arg", "family", c$message)
                       c$call <- NULL
                       stop(c)
                     })
  link <- if(!is.null(link)){
    tryCatch(
      if(family == "binomial"){
        match.arg(link, c("logit", "probit", "cauchit", "log", "cloglog"))
        } else if(family == "gaussian"){
          match.arg(link, c("identity", "log", "inverse"))
          } else if(family %in% c("negbin", "poisson")){
            match.arg(link, c("log", "sqrt", "identity"))
          },
      error = function(c){
        c$message <- gsub("'arg'", paste("'link' for", family, "family"),
                          c$message)
        c$call <- NULL
        stop(c)
      })
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
  mf <- stats::model.frame(form, data = data, na.action = stats::na.omit)
  # Correct non-standard column names
  names(mf) <- make.names(names(mf))
  if(!is.null(test_data)){
    test_data <- stats::model.frame(form, data = test_data,
                                    na.action = stats::na.omit)
    names(test_data) <- make.names(names(test_data))
    if(!all(names(mf) %in% names(test_data)))
      stop("'test_data' must contain same variables as 'data'")
  }
  n_drop <- nrow(data) - nrow(mf)
  if(n_drop > 0)
    warning(paste("Dropping", n_drop, "rows with missing data."),
            immediate. = TRUE)
  response <- names(mf)[1]
  y <- mf[,1]
  if(family == "binomial") y <- factor(y)

  #==================================================================
  # Screen for linear dependencies among predictors
  #------------------------------------------------------------------
  mm <- stats::model.matrix(form, data = mf)
  colinear_vars <- caret::findLinearCombos(mm)
  if(!is.null(colinear_vars$remove)){
    mf_to_mm <- rep(1, ncol(mf))
    factor_idx <- which(sapply(mf, class) == "factor")
    if(length(factor_idx) != 0){
      factor_exp <- sapply(mf[, factor_idx], function(x) length(levels(x)))
      mf_to_mm[factor_idx] <- factor_exp
    }
    mf_to_mm <- cumsum(mf_to_mm)
    to_remove <- names(mf)[mf_to_mm %in% colinear_vars$remove]
    stop(
      if(length(to_remove) == 1){
        paste("Linear dependency found. Consider removing predictor `",
              to_remove, "`.", sep = "")
      } else {
        paste(length(to_remove), " linear dependencies found. ",
               "Consider removing the following predictors:\n\t",
               paste0(to_remove, collapse = "\n\t"),
               sep = "")
      }
    )
  }
  #==================================================================
  # Check that number of predictors and cv folds is acceptable
  #------------------------------------------------------------------
  p <- min(ncol(mf) - 1, p_max)
  if(family == "binomial"){
    n <- min(sum(y == levels(y)[1]), sum(y == levels(y)[2]))
  } else {
    n <- nrow(mf)
  }
  alt_folds <- n / (n - p * 2)
  alt_p <- p
  while(!between(alt_folds, 1, 10)){
    alt_p <- alt_p - 1
    alt_folds <- n / (n - alt_p * 2)
  }
  if(alt_p < 1){
    if(family == "binomial"){
      stop("Your sample size for the minority class is too small.")
    } else {
      stop("Your sample size is too small.")
    }
  }
  if(alt_p < p){
    p <- alt_p
    warning(paste("'p_max' argument too high for your sample size",
                  ".\n  Reducing maximum subset size to ", p, ".",
                  sep = ""), immediate. = TRUE)
  }
  if(n_folds < alt_folds){
    n_folds <- as.integer(alt_folds)
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
  pred <- unlist(sapply(pred, function(vars){
    sapply(vars, function(x) paste0(x, collapse = " + "))}))
  pred <- c("1", pred)
  form_list <- paste(response, "~", pred)

  #======================================================================
  # Determine number of predictors for each model in list
  #----------------------------------------------------------------------
  n_pred <- sapply(pred[-1], function(x){
    x <- unlist(strsplit(x, split = " + ", fixed = TRUE))
    sum(sapply(x, function(y){
      if(is.factor(mf[,y]))
        length(levels(mf[,y])) - 1
      else
        1
    }))
  })
  n_pred <- c(intercept = 0, n_pred)

  #======================================================================
  # Eliminate formulae that exceed p_max
  #----------------------------------------------------------------------
  form_list <- form_list[n_pred <= p]
  n_pred <- n_pred[n_pred <= p]

  #======================================================================
  # Obtain fit statistics for every model
  #----------------------------------------------------------------------
  cl <- parallel::makeCluster(n_cores)
  parallel::clusterExport(cl, c("mf", "family", "dist", "link", "test_data",
                                "glm_nb", "predict_metrics",
                                "predict_metrics_", "fit_glm"),
                          envir=environment())
  fit_stats <- parallel::parLapplyLB(cl, form_list, function(form){
    fit <- fit_glm(mf, form, family, link)
    list(aic = stats::AIC(fit),
         dev = stats::deviance(fit),
         mae = mean(abs(stats::residuals(fit, type = "response"))),
         mce = -stats::logLik(fit)/nrow(mf),
         mse = mean(stats::residuals(fit, type = "response")^2),
         r2 = 1 - fit$deviance / fit$null.deviance)
    })
  stat_names <- names(fit_stats[[1]])
  get_stat <- function(stat_name){
    vapply(fit_stats, function(x) x[[stat_name]], 0.0)
  }
  fit_stats <- lapply(stat_names, get_stat)
  names(fit_stats) <- stat_names
  fit_stats <- dplyr::bind_cols(
    n_pred = n_pred,
    form = form_list,
    tibble::as_data_frame(fit_stats)
    )

  #======================================================================
  # Obtain model with maximum likelihood for each number of parameters
  #----------------------------------------------------------------------
  cv_stats <- group_by(fit_stats, n_pred) %>%
    filter(mce == min(mce)) %>%
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
      fit <- if(family == "negbin"){
        glm_nb(form, mf[i,], link = link, ...)
        } else {
          stats::glm(form, do.call(family, list(link = link)), mf[i,], ...)
        }
      stats <- try(predict_metrics(fit, test_data = mf[-i,]), silent = TRUE)
      if(class(stats) == "prediction_metrics") stats
      else list(deviance = NA_real_,
                mean_absolute_error = NA_real_,
                mean_cross_entropy = NA_real_,
                mean_squared_error = NA_real_,
                R_squared = NA_real_)
    })
  }, form_list = cv_stats$form)
  parallel::stopCluster(cl)

  #======================================================================
  # Derive cross-validation statistics
  #----------------------------------------------------------------------
  metrics <- lapply(seq_along(cv_stats$n_pred), function(i)
    transpose(at_depth(metrics, 1, i)))
  cv_mean <- at_depth(metrics, 2, function(x)
    mean(as_vector(x), na.rm = TRUE)) %>%
    transpose() %>%
    at_depth(1, as_vector) %>%
    as.data.frame()
  names(cv_mean) <- c("dev", "mae", "mce", "mse", "r2")
  cv_se <- at_depth(metrics, 2, function(x){
    x <- as_vector(x)
    sd(x, na.rm = TRUE) / sqrt(length(x[!is.na(x)]))
  }) %>% transpose() %>%
    at_depth(1, as_vector) %>%
    as.data.frame()
  names(cv_se) <- c("dev_SE", "mae_SE", "mce_SE", "mse_SE", "r2_SE")
  cv_stats <- bind_cols(cv_stats, cv_mean, cv_se)

  #======================================================================
  # Compute prediction statistics for independent test set
  #----------------------------------------------------------------------
  test_stats <- NULL
  if(!is.null(test_data)){
    metrics <- lapply(cv_stats$form, function(form){
      fit <- if(family == "negbin")
        glm_nb(form, mf, link = link, ...)
      else
        stats::glm(form, do.call(family, list(link = link)), mf, ...)
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
  structure(list(best_aic = best_aic,
                 cv_params = list(n_folds = n_folds,
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

