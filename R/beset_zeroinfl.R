#' Best Subset Selection for Zero-inflated Count Data Regression
#'
#' \code{beset_zeroinfl} performs best subset selection using repeated
#' cross-validation to find the optimal number of predictors for zero-inflated
#' regression models.
#'
#' \code{beset_zeroinfl} performs best subset selection for zero-inflated
#' count models, fitting both a count model (e.g., negative binomial regression)
#' and zero-inflation model (e.g., logistic regression) for every possible
#' combination of predictors. For models with the same number of parameters,
#' e.g., \code{m} predictors of the count process + \code{n} predictors of the
#' zero-inflation process, \code{beset_zeroinfl} picks the model with the
#' maximum log-likelihood. This results in a best fit for every \code{m * n}
#' possible number of predictors. \code{beset_zeroinfl} then uses \code{k}-fold
#' cross-validation to select the "best of the best": the best model with the
#' \code{m * n} number of predictors that minimizes prediction error, i.e., how
#' well the best models trained using \eqn{k - 1} folds predict the out-of-fold
#' sample.
#'
#' \code{beset_zeroinfl} uses \code{\link[caret]{createFolds}} to randomly
#' assign observations to \code{k} folds within subgroups based on percentiles
#' of the outcome is numeric. This insures that the folds will be matched in
#' terms of the outcome's frequency distribution. \code{beset_zeroinfl} also
#' insures the reproducibility of your analysis by requiring a \code{seed} to
#' the random number generator as one of its arguments.
#'
#' @section Warnings:
#' \enumerate{
#'  \item \code{beset_zeroinfl} handles missing data by performing listwise
#'    deletion. No other options for handling missing data are provided at this
#'    time. The user is encouraged to deal with missing values prior to running
#'    this function.
#'  \item \code{beset_zeroinfl} is intended for use with additive models only.
#'    While there is no prohibition against the inclusion of interaction or
#'    polynomial terms, this practice is strongly discouraged. At best, this
#'    will result in an inefficient search because \code{beset_zeroinfl}
#'    performs an exhaustive search over all possible variable subsets,
#'    including subsets that are hierarchically incomplete, i.e. subsets that
#'    contain an interaction term but are missing one or more of the subterms
#'    that comprise it. At worst, it may return one of these hierarchically
#'    incomplete models as the best model, an undesirable result if one cares
#'    about interpretability. If one wishes the model search to include
#'    interaction and/or non-linear effects, the
#'    \href{https://cran.r-project.org/web/packages/earth/index.html}{MARS}
#'    technique is recommended instead. \item \code{beset_zeroinfl} is best
#'    suited for searching over a small number of predictors (less than 10). For
#'    a large number of predictors (more than 20), \code{\link{beset_elnet}} is
#'    recommended instead. However, note that \code{\link{beset_elnet}} only
#'    works with a more restricted set of distributions.
#' }
#'
#' @seealso \code{\link[caret]{createFolds}}, \code{\link{beset_lm}},
#' \code{\link[pscl]{zeroinfl}}, \code{\link[base]{set.seed}},
#' \code{\link{cross_entropy}}
#'
#' @param family Character string naming the count model family. Options are
#' \code{"poisson"} (default), \code{"negbin"}, or \code{"geometric"}. (Note a
#' log link is always used with the count model).
#'
#' @param link Character string naming the link function in the binary
#' zero-inflation model. Options are \code{"logit"} (default), \code{"probit"},
#' \code{"cloglog"}, \code{"cauchit"}, or \code{"log"}. (Note a binomial family
#' is always used with the zero-inflation model).
#'
#' @param ... Optional named arguments accepted by \code{\link[pscl]{zeroinfl}}.
#'
#' @inheritParams beset_glm
#'
#' @return A "beset_zeroinfl" object with the following components:
#' \enumerate{
#'  \item\describe{
#'    \item{best_model}{an object of class \code{\link[pscl]{zeroinfl}}
#'    corresponding to the best model with the number of parameters with the
#'    smallest cross-validation error (largest cross-validation R-squared)}
#'    }
#'  \item\describe{
#'    \item{best_model_1SE}{an object of class \code{\link[pscl]{zeroinfl}}
#'    corresponding to the best model with the smallest number of
#'    parameters within one standard error of the smallest cross-validation
#'    error (largest cross-validation R-squared)}
#'    }
#'  \item\describe{
#'    \item{all_count_subsets}{a data frame containing fit statistics for every
#'      possible combination of predictors:
#'      \describe{
#'      \item{n_count_pred}{the number of predictors in model}
#'      \item{form}{formula for model}
#'      \item{train_CE}{Proportion of variance or deviance in the
#'        \code{train_data} explained by each size of best model}
#'      \item{test_CE}{if \code{test_data} is provided, the R-squared when
#'        the model fit to \code{train_data} is applied to the \code{test_data}}
#'       }
#'    }
#'  }
#'   \item\describe{
#'    \item{all_zero_subsets}{a data frame containing fit statistics for every
#'      possible combination of predictors:
#'      \describe{
#'      \item{n_zero_pred}{the number of predictors in model}
#'      \item{form}{formula for model}
#'      \item{train_CE}{Proportion of variance or deviance in the
#'        \code{train_data} explained by each size of best model}
#'      \item{test_CE}{if \code{test_data} is provided, the R-squared when
#'        the model fit to \code{train_data} is applied to the \code{test_data}}
#'       }
#'    }
#'  }
#'  \item\describe{
#'    \item{best_subsets}{a data frame containing cross-validation statistics
#'      for the best model for each combination of \code{n_count_pred} and
#'        \code{n_zero_pred}:
#'      \describe{
#'      \item{n_count_pred}{the number of predictors in model}
#'      \item{n_zero_pred}{the number of predictors in model}
#'      \item{form}{formula for best model of \code{n_count_pred} and
#'        \code{n_zero_pred}}
#'      \item{train_CE}{Proportion of variance or deviance in the
#'        \code{train_data} explained by each size of best model}
#'      \item{test_CE}{if \code{test_data} is provided, the R-squared when
#'        the model fit to \code{train_data} is applied to the \code{test_data}}
#'      \item{cv_CE}{the mean cross-validation R-squared for each size of best
#'        model, i.e., on average, how well models fit to \code{n-1} folds
#'        explain the left-out fold}
#'      \item{cv_CE_SE}{the standard error of the cross-validation R-squared for
#'       each size of best model}
#'       }
#'      }
#'    }
#'  }
#'
#' @export

beset_zeroinfl <- function(form, train_data, test_data = NULL,
                           family = "poisson", link = "logit", ...,
                           p_count_max = 10, p_zero_max = 10,
                           n_folds = 10, n_repeats = 10,
                           n_cores = NULL, seed = 42){
  #==================================================================
  # Check family and link arguments
  #------------------------------------------------------------------
  family <- try(match.arg(family, c("negbin", "poisson", "geometric")),
                silent = TRUE)
  if(class(family) == "try-error") stop("Invalid 'family' argument.")
  if(!is.null(link)){
      link <- try(match.arg(link, c("logit", "probit", "cauchit",
                                    "log", "cloglog")), silent = TRUE)
    }
  if(class(link) == "try-error") stop("Invalid 'link' argument.")

  #==================================================================
  # Screen for linear dependencies among predictors
  #------------------------------------------------------------------
  mm <- model.matrix(form, data = train_data)
  colinear_vars <- caret::findLinearCombos(mm[, 2:ncol(mm)])
  if(!is.null(colinear_vars$remove)){
    factor_idx <- which(sapply(mf, class) == "factor")
    factor_exp <- sapply(mf[, factor_idx], function(x) length(levels(x))) - 1
    mf_to_mm <- rep(1, ncol(mf))
    mf_to_mm[factor_idx] <- factor_exp
    mf_to_mm <- cumsum(mf_to_mm) - 1
    to_remove <- names(mf)[mf_to_mm %in% colinear_vars$remove]
    stop(paste(length(to_remove), " linear dependencies found. ",
               "Consider removing the following predictors:\n\t",
               paste0(to_remove, collapse = "\n\t"),
               sep = ""))
  }

  #==================================================================
  # Create model frame and extract response name and vector
  #------------------------------------------------------------------
  mf <- model.frame(form, data = train_data, na.action = na.omit)
  n_drop <- nrow(train_data) - nrow(mf)
  if(n_drop > 0)
    warning(paste("Dropping", n_drop, "rows with missing data."),
            immediate. = TRUE)
  response <- names(mf)[1]
  y <- mf[,1]
  if(min(y) != 0)
    stop("Observed lower bound does not equal 0.")

  #==================================================================
  # Check that number of predictors and cv folds is acceptable
  #------------------------------------------------------------------
  if(ncol(mf) > 11)
    stop("This approach not recommended for use with more than 10 predictors.")
  p_count <- min(ncol(mf) - 1, p_count_max)
  n <- nrow(mf)
  alt_folds <- n / (n - p_count * 10)
  alt_p_count <- p_count
  while(!dplyr::between(alt_folds, 1, 10)){
    alt_p_count <- alt_p_count - 1
    alt_folds <- n / (n - alt_p_count * 10)
  }
  if(alt_p_count < 1){
    stop("Sample size is too small.")
  }
  if(alt_p_count < p_count){
    p_count <- alt_p_count
    warning(paste("'p_count_max' argument too high given sample size",
                  ".\n Reducing maximum subset size for count-model to ",
                  p_count, ".", sep = ""), immediate. = TRUE)
  }
  if(n_folds < alt_folds){
    n_folds <- as.integer(alt_folds)
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
  alt_folds <- n / (n - p_zero * 10)
  alt_p_zero <- p_zero
  while(!dplyr::between(alt_folds, 1, 10)){
    alt_p_zero <- alt_p_zero - 1
    alt_folds <- n / (n - alt_p_zero * 10)
  }
  if(alt_p_zero < 1){
      stop("Not enough 0 values to model zero-inflation process.")
  }
  if(alt_p_zero < p_zero){
    p_zero <- alt_p_zero
    warning(paste("'p_zero_max' argument too high given counts of zero values",
                  ".\n Reducing maximum subset size for zero-model to ",
                  p_zero, ".", sep = ""), immediate. = TRUE)
  }
  if(n_folds < alt_folds){
    n_folds <- as.integer(alt_folds)
    warning(paste("'n_folds' argument too low given counts of zero values",
                  "and choice of 'p_zero_max'",
                  ".\n  Increasing number of cv folds to ", n_folds, ".",
                  sep = ""), immediate. = TRUE)
  }

  #======================================================================
  # Make list of all possible formulas with number of predictors <= p_max
  #----------------------------------------------------------------------
  n_count_pred <- unlist(sapply(1:p_count,
                                function(i) rep(i, choose(ncol(mf)-1, i))))
  n_count_pred <- c(0, n_count_pred)
  count_pred <- lapply(1:p_count, function(x)
    combn(names(mf)[2:ncol(mf)], x, simplify = FALSE))
  count_pred <- unlist(sapply(count_pred, function(vars){
    sapply(vars, function(x) paste0(x, collapse = " + "))}))
  count_pred <- c("1", count_pred)
  count_form_list <- paste(response, "~", count_pred, "| 1")
  n_zero_pred <- unlist(sapply(1:p_zero,
                                function(i) rep(i, choose(ncol(mf)-1, i))))
  n_zero_pred <- c(0, n_zero_pred)
  zero_pred <- lapply(1:p_zero, function(x)
    combn(names(mf)[2:ncol(mf)], x, simplify = FALSE))
  zero_pred <- unlist(sapply(zero_pred, function(vars){
    sapply(vars, function(x) paste0(x, collapse = " + "))}))
  zero_pred <- c("1", zero_pred)
  zero_form_list <- paste(response, "~", zero_pred, "| 1")

  #==================================================================
  # Obtain best subsets for zero vs. non-zero predictors
  #------------------------------------------------------------------
  if(is.null(n_cores)) n_cores <- parallel::detectCores() %/% 2
  cl <- parallel::makeCluster(n_cores)
  parallel::clusterExport(cl, c("mf", "family", "link", "test_data",
                                "cross_entropy"),
                          envir=environment())
  catch <- parallel::clusterEvalQ(cl, library(pscl))
  CE <- parallel::parSapplyLB(cl, zero_form_list, function(form){
    fit <- try(suppressWarnings(zeroinfl(formula(form),
                                         data = mf,
                                         dist = family,
                                         link = link,
                                         ...)), silent = TRUE)
    train_CE <- NA_real_
    test_CE <- NA_real_
    if(class(fit) == "zeroinfl"){
      train_CE <- cross_entropy(fit, mf)
      if(!is.null(test_data)) test_CE <- cross_entropy(fit, test_data)
    }
    c(train_CE, test_CE)
  })
  all_zero_subsets <- dplyr::data_frame(n_zero_pred = n_zero_pred,
                                   zero_pred = zero_pred,
                                   train_CE = CE[1,],
                                   test_CE = CE[2,])
  best_zero_subsets <- dplyr::group_by(all_zero_subsets, n_zero_pred)
  best_zero_subsets <- dplyr::filter(best_zero_subsets,
                                     train_CE == min(train_CE, na.rm = TRUE))

  #==================================================================
  # Obtain best subsets for count predictors
  #------------------------------------------------------------------
  CE <- parallel::parSapplyLB(cl, count_form_list, function(form){
    fit <- try(suppressWarnings(zeroinfl(formula(form),
                                         data = mf,
                                         dist = family,
                                         link = link,
                                         ...)), silent = TRUE)
    train_CE <- NA_real_
    test_CE <- NA_real_
    if(class(fit) == "zeroinfl"){
      train_CE <- cross_entropy(fit, mf)
      if(!is.null(test_data)) test_CE <- cross_entropy(fit, test_data)
    }
    c(train_CE, test_CE)
  })
  all_count_subsets <- dplyr::data_frame(n_count_pred = n_count_pred,
                                   count_pred = count_pred,
                                   train_CE = CE[1,],
                                   test_CE = CE[2,])
  best_count_subsets <- dplyr::group_by(all_count_subsets, n_count_pred)
  best_count_subsets <- dplyr::filter(best_count_subsets,
                                      train_CE == min(train_CE, na.rm = T))

  #======================================================================
  # Make formulas for all combinations of best zero and best count pred
  #----------------------------------------------------------------------
  best_subsets <- cbind(
    expand.grid(n_count_pred = best_count_subsets$n_count_pred,
                n_zero_pred = best_zero_subsets$n_zero_pred),
    expand.grid(count_pred = best_count_subsets$count_pred,
                zero_pred = best_zero_subsets$zero_pred,
                stringsAsFactors = FALSE),
    stringsAsFactors = FALSE)
  best_subsets$form <- paste(response, "~", best_subsets$count_pred,
                             "|", best_subsets$zero_pred)
  best_subsets <- dplyr::select(best_subsets, -count_pred, -zero_pred)


  CE <- parallel::parSapplyLB(cl, best_subsets$form, function(form){
    fit <- try(suppressWarnings(zeroinfl(formula(form),
                                         data = mf,
                                         dist = family,
                                         link = link, ...)), silent = TRUE)
    train_CE <- NA_real_
    test_CE <- NA_real_
    if(class(fit) == "zeroinfl"){
      train_CE <- cross_entropy(fit, mf)
      if(!is.null(test_data)) test_CE <- cross_entropy(fit, test_data)
    }
    c(train_CE, test_CE)
  })
  parallel::stopCluster(cl)
  best_subsets$train_CE <- CE[1,]
  best_subsets$test_CE <- CE[2,]

  #======================================================================
  # Perform cross-validation on best models
  #----------------------------------------------------------------------
  search_grid <- expand.grid(fold = 1:n_folds,
                             form = best_subsets$form,
                             stringsAsFactors = FALSE)
  seed_seq <- seq.int(from = seed, length.out = n_repeats)
  doParallel::registerDoParallel()
  CE <- foreach(seed = seed_seq, .combine = rbind, .packages = "pscl") %dopar% {
    set.seed(seed)
    folds <- vector("integer", length(y))
    # make separate fold assignments for zero vs. non-zero values to insure
    # both 0 and count process is represented in every fold
    sample_folds <- rep(1:n_folds, ceiling(length(y)/n_folds))
    folds[y == 0] <- sample(sample_folds, sum(y == 0))
    folds[y > 0] <- caret::createFolds(y[y > 0], k = n_folds, list = FALSE)
    fits <- mapply(function(fold, form){
      try(suppressWarnings(zeroinfl(formula(form),
          data = mf[folds != fold,], dist = family, link = link, ...)),
          silent = TRUE)
      }, fold = search_grid$fold, form = search_grid$form, SIMPLIFY = FALSE)

    .cross_entropy <- function(object, data){
      if(class(object) == "zeroinfl") cross_entropy(object, data) else NA_real_
    }

    cv_CE <- matrix(mapply(function(fit, fold)
        .cross_entropy(fit, mf[folds == fold,]),
        fit = fits, fold = search_grid$fold),
        nrow = n_folds, ncol = nrow(search_grid))
    CE_cv <- apply(cv_CE, 2, mean, na.rm = T)
    CE_cv_SE <- apply(cv_CE, 2, function(x) sqrt(var(x, na.rm = T)/length(x)))

    data.frame(form = search_grid$form,
               cv_CE = CE_cv,
               cv_CE_SE = CE_cv_SE,
               stringsAsFactors = FALSE)
  }
  #======================================================================
  # Derive cross-validation statistics
  #----------------------------------------------------------------------
  CE <- dplyr::group_by(CE, form)
  mean_CE <- dplyr::summarize_each(CE, dplyr::funs(mean(., na.rm = T)))
  best_subsets <- dplyr::left_join(best_subsets, mean_CE, by = "form")

  #======================================================================
  # Determine best model and best model within 1SE with smallest number of
  # total predictors (number of zero predictors + number of count predictors)
  #----------------------------------------------------------------------

  best <- which.min(best_subsets$cv_CE)
  best_form <- best_subsets$form[best]
  max_CE <- min(best_subsets$cv_CE) + best_subsets$cv_CE_SE[best]
  best_subsets_1SE <- best_subsets[best_subsets$cv_CE < max_CE,]
  best_subsets_1SE$total_pred <- with(best_subsets_1SE,
                                      n_count_pred + n_zero_pred)
  best_subsets_1SE <- dplyr::group_by(best_subsets_1SE, total_pred)
  best_subsets_1SE <- dplyr::filter(best_subsets_1SE, cv_CE == min(cv_CE))
  best_form_1SE <- best_subsets_1SE$form[which.min(best_subsets_1SE$total_pred)]

  best_model <- pscl::zeroinfl(formula(best_form), mf, dist = family,
                               link = link, ...)
    if(best_form == best_form_1SE){
      best_model_1SE <- best_model
    } else {
      best_model_1SE <- pscl::zeroinfl(formula(best_form_1SE), mf,
                                       dist = family, link = link, ...)
    }

  #======================================================================
  # Construct beset_zeroinfl object
  #----------------------------------------------------------------------

  structure(list(all_zero_subsets = all_zero_subsets,
                 all_count_subsets = all_count_subsets,
                 best_subsets = best_subsets,
                 best_model = best_model, best_model_1SE = best_model_1SE),
          class = "beset_zeroinfl")
}
