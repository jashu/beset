#' Best Subset Selection for Generalized Linear Models
#'
#' \code{beset_glm} performs best subset selection using repeated
#' cross-validation to find the optimal number of predictors.
#'
#' \code{beset_glm} performs best subset selection for generalized linear
#' models, fitting a separate model for each possible combination of predictors,
#' i.e., all models that contain exactly 1 predictor, all models that contain
#' exactly 2 predictors, and so forth. For each number of predictors,
#' \code{beset_glm} picks the model with the smallest deviance. This results in
#' a best fit for every possible number of predictors. \code{beset_glm} then
#' uses \code{k}-fold cross-validation to select the "best of the best": the
#' best model with the number of predictors that minimizes prediction error,
#' i.e., how well the best models trained using \eqn{k - 1} folds predict the
#' out-of-fold sample. \code{beset_glm} uses \code{\link[caret]{createFolds}} to
#' randomly assign observations to \code{k} folds within levels of the outcome
#' when the outcome is a factor or within subgroups based on percentiles when
#' the outcome is numeric. This insures that the folds will be matched in terms
#' of the outcome's frequency distribution. \code{beset_glm} also insures the
#' reproducibility of your analysis by requiring a \code{seed} to the random
#' number generator as one of its arguments.
#'
#' @section Warnings:
#' \enumerate{
#'  \item \code{beset_glm} handles missing data by performing listwise deletion.
#'   No other options for handling missing data are provided. The user is
#'   encouraged to deal with missing values prior to running this function.
#'  \item \code{beset_glm} is intended for use with additive models only.
#'    While there is no prohibition against the inclusion of interaction or
#'    polynomial terms, this practice is strongly discouraged. At best, this
#'    will result in an inefficient search because \code{beset_glm} performs an
#'    exhaustive search over all possible variable subsets, including subsets
#'    that are hierarchically incomplete, i.e. subsets that contain an
#'    interaction term but are missing one or more of the subterms that comprise
#'    it. At worst, it may return one of these hierarchically incomplete models
#'    as the best model, an undesirable result if one cares about
#'    interpretability. If one wishes the model search to include interaction
#'    and/or non-linear effects, the
#'    \href{https://cran.r-project.org/web/packages/earth/index.html}{MARS}
#'    technique is recommended instead. \item \code{beset_glm} is best suited
#'    for searching over a small number of predictors (less than 10). For a
#'    large number of predictors (more than 20), \code{\link{beset_elnet}} is
#'    recommended instead. However, note that \code{\link{beset_elnet}} only
#'    works with a more restricted set of distributions.
#' }
#'
#' @seealso \code{\link[caret]{createFolds}}, \code{\link{beset_lm}},
#' \code{\link[base]{set.seed}}, \code{\link{predict_R2}}
#'
#' @param form A model \code{\link[stats]{formula}}.
#'
#' @param train_data A \code{\link[base]{data.frame}} with the variables in
#' \code{form} and the data to be used in model training.
#'
#' @param test_data Optional \code{\link[base]{data.frame}} with the variables
#' in \code{form} and the data to be used in model testing.
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
#' @return A list with the following three components:
#' \enumerate{
#'  \item\describe{
#'    \item{best_model}{an object of class \code{\link[stats]{glm}} corresponding
#'    to the best model with the number of parameters with the smallest
#'    cross-validation error (largest cross-validation R-squared)}
#'    }
#'  \item\describe{
#'    \item{best_model_1SE}{an object of class \code{\link[stats]{glm}}
#'    corresponding to the best model with the smallest number of
#'    parameters within one standard error of the smallest cross-validation
#'    error (largest cross-validation R-squared)}
#'    }
#'  \item \code{\link[base]{data.frame}} with the following columns:
#'    \describe{
#'    \item{k_folds}{the number of k-folds used}
#'    \item{n_preds}{the number of predictors in model}
#'    \item{R2_train}{the mean in-fold (training) R-squared for each size of
#'      best model}
#'    \item{R2_train_SE}{the standard error of the mean in-fold (training)
#'      R-squared for each size of best model}
#'    \item{R2_cv}{the mean out-of-fold (test) R-squared for each size of best
#'      model}
#'    \item{R2_cv_SE}{the standard error of the mean out-of-fold (test)
#'      R-squared for each size of best model}
#'    \item{R2_test}{if \code{test_data} is provided, the mean R-squared when
#'      the model fit to each subsample of the training data is applied to the
#'      independent test set}
#'    \item{R2_test_SE}{if \code{test_data} is provided, the standard error of
#'      the mean R-squared when the model fit to each subsample of the training
#'      data is applied to the independent test set}
#'    }
#'  }
#'
#' @export

beset_glm <- function(form, train_data, test_data = NULL,
                      family = "gaussian", ..., p_max = 10,
                      n_folds = 10, n_repeats = 10,
                      n_cores = NULL, seed = 42){
  if(family == "gaussian") return(beset_lm(form, train_data, test_data, p_max,
                                           n_folds, n_repeats, seed))
  mf <- model.frame(form, data = train_data)
  n_drop <- nrow(train_data) - nrow(mf)
  if(n_drop > 0) warning(paste("Dropping", n_drop, "rows with missing data."),
                         immediate. = TRUE)
  p <- min(ncol(mf) - 1, p_max)
  if(p > 20){
    stop(paste("`beset_glm` does not allow more than 20 predictors."))
  }
  n <- nrow(mf)
  alt_p <- n - 1 - n %/% n_folds
  if(p > alt_p){
    warning(paste("You have more predictors than your sample size will support",
                  ".\n  Setting maximum subset size to ", alt_p, ".",
                  sep = ""), immediate. = TRUE)
    p <- alt_p
  }
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
  get_se <- function(x) sqrt(var(x)/length(x))
  y <- mf[,1]
  response <- names(mf)[1]
  R2 <- expand.grid(n_repeats = 1:n_repeats,
                    n_preds = 1:p,
                    R2_train = NA_real_,
                    R2_train_SE = NA_real_,
                    R2_cv = NA_real_,
                    R2_cv_SE = NA_real_,
                    R2_test = NA_real_,
                    R2_test_SE = NA_real_)
  all_subsets <- as.list(1:p)
  all_subsets <- lapply(all_subsets, function(x)
    combn(2:ncol(mf), x, simplify = FALSE))
  search_list <- as.list(1:(n_folds * p))
  search_grid <- expand.grid(folds = 1:n_folds, n_pred = 1:p)
  search_list <- lapply(search_list, function(x) as.list(search_grid[x,]))
  pb <- txtProgressBar(min = 0, max = n_repeats, style = 3)
  counter <- 1
  set.seed(seed)
  for(n in 1:n_repeats){
    folds <- caret::createFolds(y, k = n_folds, list = FALSE)
    all_fits <- lapply(search_list, function(search_param){
      lapply(all_subsets[[search_param$n_pred]], function(col_id){
        glm(paste(response, ".", sep = "~"),
            data = mf[folds != search_param$folds, c(1, col_id)],
            family = family)
        })
      })
    best_fits <- lapply(all_fits, function(x){
      x[[which.min(sapply(x, deviance))]]
    })

    idx <- R2$n_repeats == n

    train_R2 <- matrix(sapply(best_fits, function(x)
      1 - x$deviance / x$null.deviance), nrow = n_folds, ncol = p)
    R2$R2_train[idx] <- apply(train_R2, 2, mean)
    R2$R2_train_SE[idx] <- apply(train_R2, 2, get_se)

    cv_R2 <- matrix(mapply(function(fit, fold)
      predict_R2(fit, mf[folds == fold,]),
      fit = best_fits, fold = search_grid$folds), nrow = n_folds, ncol = p)
    R2$R2_cv[idx] <- apply(cv_R2, 2, mean)
    R2$R2_cv_SE[idx] <- apply(cv_R2, 2, get_se)

    if(!is.null(test_data)){
      test_R2 <- matrix(sapply(best_fits, function(x)
        predict_R2(x, test_data)), nrow = n_folds, ncol = p)
      R2$R2_test[idx] <- apply(test_R2, 2, mean)
      R2$R2_test_SE[idx] <- apply(test_R2, 2, get_se)
    }

    setTxtProgressBar(pb, counter)
    counter <- counter + 1
  }
  R2 <- dplyr::group_by(R2, n_preds)
  R2 <- dplyr::summarize_each(R2, dplyr::funs(mean(., na.rm = T)))
  max_R2 <- max(R2$R2_cv)
  if(max_R2 <= 0){
    best_model <- glm(paste(response, "~ 1"), data = mf, family = family, ...)
    best_model_1SE <- best_model
  } else {
    n_pred <- R2$n_preds[R2$R2_cv == max_R2]
    all_fits <- lapply(all_subsets[[n_pred]],
                function(col_idx) glm(paste(response, ".", sep = "~"),
                                      data = mf[, c(1, col_idx)],
                                      family = family, ...))
    best_model <- all_fits[[which.min(sapply(all_fits, deviance))]]
    min_R2 <- max_R2 - R2$R2_cv_SE[R2$R2_cv == max_R2]
    if(min_R2 <= 0){
      best_model_1SE <- glm(paste(response, "~ 1"), data = mf,
                            family = family, ...)
    } else {
      temp <- R2[R2$R2_cv > min_R2,]
      n_pred <- min(temp$n_preds)
      all_fits <- lapply(all_subsets[[n_pred]],
                         function(col_idx) glm(paste(response, ".", sep = "~"),
                                               data = mf[, c(1, col_idx)],
                                               family = family, ...))
      best_model_1SE <- all_fits[[which.min(sapply(all_fits, deviance))]]
    }
  }
  structure(list(R2 = R2, best_model = best_model,
                 best_model_1SE = best_model_1SE),
            class = "beset_glm")
}
