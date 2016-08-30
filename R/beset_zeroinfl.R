#' Best Subset Selection for Generalized Linear Models
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
#' \code{\link{predict_R2}}
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
#' @param search Character string naming the search strategy to use. Options are
#' \code{"stagewise"} (default), \code{"identical"}, or \code{"exhasutive"}. See
#' 'Details' for explanation.
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
#'    \item{all_subsets}{a data frame containing fit statistics for every
#'      possible combination of predictors:
#'      \describe{
#'      \item{n_pred}{the number of predictors in model}
#'      \item{form}{formula for model}
#'      \item{train_R2}{Proportion of variance or deviance in the
#'        \code{train_data} explained by each size of best model}
#'      \item{test_R2}{if \code{test_data} is provided, the R-squared when
#'        the model fit to \code{train_data} is applied to the \code{test_data}}
#'       }
#'    }
#'  }
#'  \item\describe{
#'    \item{best_subsets}{a data frame containing cross-validation statistics
#'      for the best model for each \code{n_pred} listed in \code{all_subsets}:
#'      \describe{
#'      \item{n_pred}{the number of predictors in model}
#'      \item{form}{formula for best model of \code{n_pred}}
#'      \item{train_R2}{Proportion of variance or deviance in the
#'        \code{train_data} explained by each size of best model}
#'      \item{test_R2}{if \code{test_data} is provided, the R-squared when
#'        the model fit to \code{train_data} is applied to the \code{test_data}}
#'      \item{cv_R2}{the mean cross-validation R-squared for each size of best
#'        model, i.e., on average, how well models fit to \code{n-1} folds
#'        explain the left-out fold}
#'      \item{cv_R2_SE}{the standard error of the cross-validation R-squared for
#'       each size of best model}
#'       }
#'    }
#'  }
#'  item\describe{
#'    \item{best_binomial}{an object of class \code{\link{beset_glm}}
#'    corresponding to the best subsets of binomial models predicting zero vs.
#'    non-zero values}
#'    }
#'
#'  item\describe{
#'    \item{search}{record of the strategy used to constrain model search}
#'    }
#' }
#'
#' @export

beset_zeroinfl <- function(form, train_data, test_data = NULL,
                           family = "poisson", link = "logit", ...,
                           search = "stagewise",
                           p_max = 10, n_folds = 10, n_repeats = 10,
                           n_cores = NULL, seed = 42){
  #==================================================================
  # Check family, link, and search arguments
  #------------------------------------------------------------------
  family <- try(match.arg(family, c("negbin", "poisson", "geometric")),
                silent = TRUE)
  if(class(family) == "try-error") stop("Invalid 'family' argument.")
  if(!is.null(link)){
      link <- try(match.arg(link, c("logit", "probit", "cauchit",
                                    "log", "cloglog")), silent = TRUE)
    }
  if(class(link) == "try-error") stop("Invalid 'link' argument.")
  search <- try(match.arg(family, c("stagewise", "identical", "exhaustive")),
                silent = TRUE)
  if(class(family) == "try-error") stop("Invalid 'search' argument.")

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
  p <- min(ncol(mf) - 1, p_max)
  mf_binom <- mf
  mf_binom[,1] <- mf[,1] > 0
  y_binom <- as.factor(mf_binom[,1])
  n <- min(sum(y_binom == levels(y_binom)[1]),
           sum(y_binom == levels(y_binom)[2]))
  alt_folds <- n / (n - p * 10)
  alt_p <- p
  while(!dplyr::between(alt_folds, 1, 10)){
    alt_p <- alt_p - 1
    alt_folds <- n / (n - alt_p * 10)
  }
  if(alt_p < 1){
      stop("Sample size too small or not enough 0 values.")
  }
  if(alt_p < p){
    p <- alt_p
    warning(paste("'p_max' argument too high given counts of zero values",
                  ".\n  Reducing maximum subset size to ", p, ".",
                  sep = ""), immediate. = TRUE)
  }
  if(n_folds < alt_folds){
    n_folds <- as.integer(alt_folds)
    warning(paste("'n_folds' argument too low given counts of zero values",
                  "and the choice of 'p_max'",
                  ".\n  Increasing number of cv folds to ", n_folds, ".",
                  sep = ""), immediate. = TRUE)
  }

  #==================================================================
  # Obtain best subsets for zero vs. non-zero predictors
  #------------------------------------------------------------------
  best_binomial <- beset_glm(form, mf_binom, family = "binomial",
                             link = link, p_max = p, n_folds = n_folds,
                             n_repeats = n_repeats, n_cores = n_cores,
                             seed = seed)
  n_pred <- best_binomial$all_subsets$n_pred
  pred <- lapply(1:p, function(x)
    combn(names(mf)[2:ncol(mf)], x, simplify = FALSE))
  pred <- unlist(sapply(pred, function(vars){
    sapply(vars, function(x) paste0(x, collapse = " + "))
  }))
  pred <- c("1", pred)
  if(search == "exhaustive"){
    pred_grid <- expand.grid(count_pred = pred,
                             zero_pred = pred,
                             stringsAsFactors = FALSE)
    all_subsets <- expand.grid(n_count_pred = n_pred, n_zero_pred = n_pred)
  }
  if(search == "identical"){
    pred_grid <- data.frame(count_pred = pred,
                            zero_pred = pred,
                            stringsAsFactors = FALSE)
    all_subsets <- data.frame(n_count_pred = n_pred, n_zero_pred = n_pred)
  }
  if(search == "stagewise"){
      zero_pred_1SE <- as.character(terms(best_binomial$best_model_1SE))[3]
      n_zero_pred_1SE <- length(unlist(
        strsplit(zero_pred_1SE, " + ", fixed = TRUE)))
      zero_pred <- as.character(terms(best_binomial$best_model))[3]
      n_zero_pred <- length(unlist(strsplit(zero_pred, " + ", fixed = TRUE)))
      if(n_zero_pred != n_zero_pred_1SE){
        zero_pred <- c(zero_pred_1SE, zero_pred)
        n_zero_pred <- c(n_zero_pred_1SE, n_zero_pred)
      }
      pred_grid <- expand.grid(count_pred = pred,
                               zero_pred = zero_pred,
                               stringsAsFactors = FALSE)
      all_subsets <- expand.grid(n_count_pred = n_pred,
                               n_zero_pred = n_zero_pred)
  }
  all_subsets$form <- paste(response, "~", pred_grid$count_pred, "|",
                            pred_grid$zero_pred)

  #======================================================================
  # Obtain R^2 for every model; use parallel computing if possible
  #----------------------------------------------------------------------
  if(is.null(n_cores)) n_cores <- parallel::detectCores() %/% 2
  cl <- parallel::makeCluster(n_cores)
  parallel::clusterExport(cl, c("mf", "family", "link", "test_data",
                                "predict_R2"),
                          envir=environment())
  catch <- parallel::clusterEvalQ(cl, library(pscl))
  R2 <- parallel::parSapplyLB(cl, all_subsets$form, function(form){
    fit <- try(suppressWarnings(zeroinfl(formula(form),
                                         data = mf,
                                         dist = family,
                                         link = link,
                                         ...)),
               silent = TRUE)
    train_R2 <- NA_real_
    test_R2 <- NA_real_
    if(class(fit) == "zeroinfl"){
      train_R2 <- predict_R2(fit, mf)
      if(!is.null(test_data)) test_R2 <- predict_R2(fit, test_data)
    }
    c(train_R2, test_R2)
    })
  parallel::stopCluster(cl)
  all_subsets$train_R2 <- R2[1,]
  all_subsets$test_R2 <- R2[2,]
  all_subsets <- dplyr::arrange(all_subsets, n_zero_pred, n_count_pred,
                                dplyr::desc(train_R2))

  #======================================================================
  # Obtain model with best R^2 for each number of parameters
  #----------------------------------------------------------------------
  best_subsets <- dplyr::group_by(all_subsets, n_zero_pred, n_count_pred)
  best_subsets <- dplyr::filter(best_subsets,
                                train_R2 == max(train_R2, na.rm = TRUE))

  #======================================================================
  # Perform cross-validation on best models
  #----------------------------------------------------------------------
  unique_zero_pred <- sort(unique(best_subsets$n_zero_pred))
  search_grid <- expand.grid(fold = 1:n_folds,
                             n_count_pred = 0:p,
                             n_zero_pred = unique_zero_pred)
  search_grid <- dplyr::left_join(search_grid, best_subsets,
                                  by = c("n_count_pred", "n_zero_pred"))
  seed_seq <- seq.int(from = seed, length.out = n_repeats)
  doParallel::registerDoParallel()
  R2 <- foreach(seed = seed_seq, .combine = rbind, .packages = "pscl") %dopar% {
    set.seed(seed)
    folds <- vector("integer", length(y))
    # make separate fold assignments for zero vs. non-zero values to insure
    # both 0 and count process is represented in every fold
    folds[y == 0] <- sample.int(n_folds, sum(y == 0), replace = TRUE)
    folds[y > 0] <- caret::createFolds(y[y > 0], k = n_folds, list = FALSE)
    fits <- mapply(function(fold, form){
      try(suppressWarnings(zeroinfl(formula(form),
          data = mf[folds != fold,], dist = family, link = link, ...)),
          silent = TRUE)
      }, fold = search_grid$fold, form = search_grid$form, SIMPLIFY = FALSE)

  .predict_R2 <- function(object, data){
    if(class(object) == "zeroinfl") predict_R2(object, data) else NA_real_
  }
    cv_R2 <- matrix(mapply(function(fit, fold)
      .predict_R2(fit, mf[folds == fold,]),
      fit = fits, fold = search_grid$fold),
      nrow = n_folds, ncol = nrow(search_grid))
    R2_cv <- apply(cv_R2, 2, mean, na.rm = T)
    R2_cv_SE <- apply(cv_R2, 2, function(x) sqrt(var(x, na.rm = T)/length(x)))

    data.frame(n_count_pred = search_grid$n_count_pred,
               n_zero_pred = search_grid$n_zero_pred,
               cv_R2 = R2_cv,
               cv_R2_SE = R2_cv_SE)
  }
  #======================================================================
  # Derive cross-validation statistics
  #----------------------------------------------------------------------
  R2 <- dplyr::group_by(R2, n_count_pred, n_zero_pred)
  mean_R2 <- dplyr::summarize_each(R2, dplyr::funs(mean(., na.rm = T)))
  # if there are more repeats than folds in each repeat,
  # report standard deviation of the mean R2 across repeats (empirical SEM);
  # otherwise report mean theoretical SEM based on averaging R2 across folds.
  if(n_repeats >= n_folds){
    se_R2 <- dplyr::summarize(R2, cv_R2_SE = sd(cv_R2, na.rm = T))
    mean_R2$cv_R2_SE <- se_R2$cv_R2_SE
  }

  #======================================================================
  # Fit best model and 1SE best model to full data set
  #----------------------------------------------------------------------
  best_subsets <- dplyr::left_join(best_subsets, mean_R2,
                                   by = c("n_count_pred", "n_zero_pred"))
  best_subset <- which.max(best_subsets$cv_R2)
  best_form <- best_subsets$form[best_subset]
  min_R2 <- max(best_subsets$cv_R2) - best_subsets$cv_R2_SE[best_subset]
  best_subsets_1SE <- best_subsets[best_subsets$cv_R2 > min_R2,]
  total_pred <- best_subsets_1SE$n_zero_pred + best_subsets_1SE$n_count_pred
  best_form_1SE <- best_subsets_1SE$form[which.min(total_pred)]

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

  structure(list(all_subsets = all_subsets, best_subsets = best_subsets,
               best_model = best_model, best_model_1SE = best_model_1SE,
               best_binomial = best_binomial, search = search),
          class = "beset_zeroinfl")
}
